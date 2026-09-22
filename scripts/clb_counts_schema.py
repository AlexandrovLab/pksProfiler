#!/usr/bin/env python3
"""The contract a per-sample clb count table must meet, for either method.

F06. The HMM merger accepted an unknown gene, a negative count, a missing count
silently turned into zero, a header-only sample, fractional counts, and duplicate gene
rows quietly summed. The alignment merger checked most of that, so the two halves of
the same cohort table were validated to different standards while the release notes
claimed both enforced nineteen unique ordered genes.

One definition, used by both:

  * exactly clbA through clbS, no more and no fewer
  * one row per gene -- a repeated gene is an error, not something to add up
  * counts are finite, non-negative integers. This pipeline never asks featureCounts
    or HMMER for fractional assignment, so a fraction means something else happened
  * a sample whose nineteen counts are all zero is valid and must stay in the cohort;
    an empty or truncated file is not the same thing and is an error
"""
import math
from collections import Counter

CLB_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]
CLB_GENE_SET = set(CLB_GENES)


class CountSchemaError(ValueError):
    """A per-sample count table that cannot be trusted into the cohort matrix."""


def validate_clb_counts(rows, source):
    """
    rows: iterable of (gene, value) as read from the file, values still text.
    Returns {gene: int} covering exactly CLB_GENES, in that order.

    Raises CountSchemaError naming the file and what was wrong with it.
    """
    pairs = [(str(gene).strip(), value) for gene, value in rows
             if str(gene).strip() != ""]

    if not pairs:
        raise CountSchemaError(
            f"{source}: no gene rows. A sample with no clb reads must still carry "
            f"{len(CLB_GENES)} rows of zero; an empty file is a truncated one.")

    counted = Counter(gene for gene, _ in pairs)
    duplicated = sorted(gene for gene, n in counted.items() if n > 1)
    if duplicated:
        raise CountSchemaError(
            f"{source}: duplicate gene row(s): {', '.join(duplicated)}. Two rows for "
            "one gene is a malformed file, not a count to be summed.")

    observed = {gene for gene, _ in pairs}
    unknown = sorted(observed - CLB_GENE_SET)
    if unknown:
        raise CountSchemaError(
            f"{source}: not a clb gene: {', '.join(unknown)}. Expected exactly "
            "clbA-clbS.")

    missing = [gene for gene in CLB_GENES if gene not in observed]
    if missing:
        raise CountSchemaError(
            f"{source}: missing {len(missing)} of {len(CLB_GENES)} genes: "
            f"{', '.join(missing)}. A gene with no reads is a zero, not an absent row.")

    counts = {}
    for gene, value in pairs:
        text = str(value).strip()
        if text == "":
            raise CountSchemaError(
                f"{source}: {gene} has no count. An empty field is not a zero.")
        try:
            number = float(text)
        except ValueError:
            raise CountSchemaError(f"{source}: {gene} count is not a number: {text!r}")
        if not math.isfinite(number):
            raise CountSchemaError(f"{source}: {gene} count is not finite: {text!r}")
        if number < 0:
            raise CountSchemaError(f"{source}: {gene} count is negative: {text!r}")
        if number != int(number):
            raise CountSchemaError(
                f"{source}: {gene} count is fractional: {text!r}. Neither featureCounts "
                "nor HMMER is asked for fractional assignment here.")
        counts[gene] = int(number)

    return {gene: counts[gene] for gene in CLB_GENES}
