#!/usr/bin/env python3
"""One preflight contract: check everything the run needs, before anything runs.

F20. Validation existed but was scattered and asked the wrong question. Stage
prerequisites checked that a flag was *set* -- `--pks_taxa` requires `--kraken_db` --
and never that what it pointed at was usable. `--hg38_db /typo/human.mmi` passed every
check and failed in the first mapReads task, after every sample had been extracted.
The checks also aborted one at a time, so fixing a cohort's worth of paths meant a
launch per mistake.

This runs once at launch, before any task is submitted, and reports every problem it
finds in one pass.

    preflight.py --config preflight.json [--quiet]

Exit 1 if anything is wrong. Warnings alone do not fail the run.
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

# Written by hmmpress; nhmmscan needs them beside the model.
HMM_PRESS_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")
BOWTIE2_SUFFIXES = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
CLB_GENES = [f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS"]
SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

ALIGNMENT_SUFFIXES = (".bam", ".cram", ".sam")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


class Report:
    """Problems stop the run; notes do not. Both are reported together."""

    def __init__(self):
        self.problems = []
        self.notes = []

    def problem(self, section, message):
        self.problems.append((section, message))

    def note(self, section, message):
        self.notes.append((section, message))

    def require_path(self, section, label, value, kind="file"):
        """Present, and pointing at something that is there."""
        if not value:
            self.problem(section, f"{label} is not set")
            return None
        path = Path(str(value))
        if not path.exists():
            self.problem(section, f"{label} does not exist: {path}")
            return None
        if kind == "dir" and not path.is_dir():
            self.problem(section, f"{label} is not a directory: {path}")
            return None
        if kind == "file" and not path.is_file():
            self.problem(section, f"{label} is not a file: {path}")
            return None
        return path


def check_sample_sheet(config, report):
    section = "sample sheet"
    sheet = report.require_path(section, "--sample", config.get("sample"))
    if not sheet:
        return

    with sheet.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        report.problem(section, f"contains no samples: {sheet}")
        return

    columns = set(rows[0].keys())
    input_type = str(config.get("input_data_type", "auto"))
    if "patient" not in columns:
        report.problem(section, "missing the required `patient` column")
    if input_type in ("bam", "cram", "auto") and not ({"alignment", "bam", "cram"} & columns):
        report.problem(section, f"input_data_type={input_type} needs an alignment, bam or cram column")
    if input_type == "fastq" and "fastq1" not in columns:
        report.problem(section, "input_data_type=fastq needs a fastq1 column")

    seen = {}
    for number, row in enumerate(rows, start=2):
        sample = str(row.get("patient", "")).strip()
        if not sample:
            report.problem(section, f"row {number}: empty patient value")
            continue
        if not SAMPLE_ID.match(sample):
            report.problem(section, f"row {number}: invalid patient value {sample!r}; "
                                    "letters, numbers, period, underscore and hyphen only, "
                                    "starting alphanumeric")
        if sample in seen:
            report.problem(section, f"row {number}: duplicate patient value {sample!r} "
                                    f"(already on row {seen[sample]})")
        seen[sample] = number

        for column in ("alignment", "bam", "cram", "fastq1", "fastq2"):
            value = str(row.get(column, "") or "").strip()
            if not value:
                continue
            path = Path(value)
            if not path.is_file():
                report.problem(section, f"row {number}: {column} does not exist: {value}")
                continue
            if path.stat().st_size == 0:
                report.problem(section, f"row {number}: {column} is empty: {value}")
            name = path.name.lower()
            if column in ("alignment", "bam", "cram") and not name.endswith(ALIGNMENT_SUFFIXES):
                report.note(section, f"row {number}: {column} is not .bam/.cram/.sam: {path.name}")
            if column.startswith("fastq") and not name.endswith(FASTQ_SUFFIXES):
                report.note(section, f"row {number}: {column} is not a FASTQ name: {path.name}")

        cram = str(row.get("cram", "") or "").strip()
        alignment = str(row.get("alignment", "") or "").strip()
        looks_cram = cram or alignment.lower().endswith(".cram")
        if looks_cram and not config.get("cram_reference"):
            report.problem(section, f"row {number}: CRAM input needs --cram_reference")

    report.note(section, f"{len(rows)} samples")


def check_references(config, report):
    section = "references"
    # Required for every run, and never checked before: a typo here used to surface in
    # the first mapReads task, after extraction had run for the whole cohort.
    report.require_path(section, "--hg38_db", config.get("hg38_db"))
    report.require_path(section, "--t2t_phix_db", config.get("t2t_phix_db"))
    if config.get("pangenome_db"):
        report.require_path(section, "--pangenome_db", config.get("pangenome_db"))
    if config.get("cram_reference"):
        report.require_path(section, "--cram_reference", config.get("cram_reference"))

    # Optional: unset means fastp uses its own detection. When given, a typo used to
    # surface in the first filterReads task, which is the failure F20 exists to prevent.
    if config.get("adapters"):
        report.require_path(section, "--adapters", config.get("adapters"))
    report.require_path(section, "reference FASTA", config.get("pks_reference_fasta"))
    report.require_path(section, "clb annotation", config.get("pks_genome_annotation"))

    for label, prefix in (("profiling index", config.get("pks_genome")),
                          ("targeted recruitment index", config.get("pks_recruit_index"))):
        if not prefix:
            continue
        missing = [suffix for suffix in BOWTIE2_SUFFIXES if not Path(f"{prefix}{suffix}").is_file()]
        if missing:
            report.problem(section, f"{label} {prefix} is missing "
                                    f"{len(missing)} of {len(BOWTIE2_SUFFIXES)} Bowtie2 files: "
                                    f"{', '.join(missing)}")


def check_hmm_models(config, report):
    section = "HMM models"
    if str(config.get("profiling_method", "bowtie2")) not in ("hmm", "both"):
        needs_nucleotide = False
    else:
        needs_nucleotide = True

    if needs_nucleotide:
        model = report.require_path(section, "--hmm_model", config.get("hmm_model"))
        if model:
            missing = [s for s in HMM_PRESS_SUFFIXES if not Path(f"{model}{s}").is_file()]
            if missing:
                report.problem(section, f"{model} has not been pressed: missing "
                                        f"{', '.join(missing)}. Run `hmmpress {model}`")
            present = set(re.findall(r"^NAME\s+(\S+)", model.read_text(errors="ignore"), re.M))
            absent = [gene for gene in CLB_GENES
                      if not any(gene.lower() in name.lower() for name in present)]
            if absent:
                report.problem(section, f"{model} does not cover {len(absent)} clb genes: "
                                        f"{', '.join(absent)}")

    if any(config.get(flag) for flag in ("enable_mags", "tumor_enable_mags",
                                         "tumor_full_contig_context")):
        report.require_path(section, "--clb_protein_hmm", config.get("clb_protein_hmm"))
        report.require_path(section, "clb protein FASTA", config.get("clb_protein_fasta"))


def check_taxonomy(config, report):
    section = "taxonomy"
    wants_taxonomy = any(config.get(flag) for flag in ("pks_taxa", "pks_community_taxa"))
    wants_prefilter = str(config.get("prefilter_mode", "off")) == "balanced"
    if not (wants_taxonomy or wants_prefilter):
        return

    kraken = report.require_path(section, "--kraken_db", config.get("kraken_db"), kind="dir")
    if not kraken:
        return
    for member in ("database.kdb", "database.idx"):
        if not (kraken / member).is_file():
            report.problem(section, f"--kraken_db {kraken} does not look like a KrakenUniq "
                                    f"database: {member} is missing")

    if not wants_taxonomy:
        return

    # Bracken cannot invent a distribution: the database must already carry one built
    # for this read length, and the failure otherwise comes late and unexplained.
    length = config.get("bracken_read_length")
    if not length:
        report.problem(section, "taxonomic profiling requires --bracken_read_length")
        return
    available = sorted(int(match.group(1)) for match in
                       (re.match(r"database(\d+)mers\.kmer_distrib$", p.name)
                        for p in kraken.iterdir()) if match)
    if not available:
        report.problem(section, f"--kraken_db {kraken} carries no Bracken read-length "
                                "distributions (database<N>mers.kmer_distrib)")
    elif int(length) not in available:
        report.problem(section, f"--bracken_read_length {length} has no distribution in "
                                f"{kraken}. Available: {', '.join(str(a) for a in available)}")


def check_stage_databases(config, report):
    section = "stage databases"
    if any(config.get(flag) for flag in ("enable_mags", "tumor_enable_mags")):
        report.require_path(section, "--gtdbtk_db", config.get("gtdbtk_db"), kind="dir")
        report.require_path(section, "--checkm2_db", config.get("checkm2_db"))
        report.require_path(section, "--genomad_db", config.get("genomad_db"), kind="dir")
    elif config.get("tumor_full_contig_context"):
        report.require_path(section, "--genomad_db", config.get("genomad_db"), kind="dir")

    if config.get("enable_strain_typing"):
        report.require_path(section, "strain typing lookup", config.get("st_phylogroup_lookup"))


def run_checks(config):
    report = Report()
    check_sample_sheet(config, report)
    check_references(config, report)
    check_hmm_models(config, report)
    check_taxonomy(config, report)
    check_stage_databases(config, report)
    return report


def render(report):
    lines = []
    if report.problems:
        lines.append("preflight found %d problem%s:"
                     % (len(report.problems), "" if len(report.problems) == 1 else "s"))
        for section, message in report.problems:
            lines.append(f"  [{section}] {message}")
    if report.notes:
        lines.append("notes:")
        for section, message in report.notes:
            lines.append(f"  [{section}] {message}")
    if not report.problems:
        lines.append("preflight passed")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--quiet", action="store_true", help="print only problems")
    args = parser.parse_args()

    config = json.loads(args.config.read_text())
    report = run_checks(config)
    if report.problems or not args.quiet:
        print(render(report))
    return 1 if report.problems else 0


if __name__ == "__main__":
    sys.exit(main())
