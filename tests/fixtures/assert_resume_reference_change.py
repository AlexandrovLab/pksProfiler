#!/usr/bin/env python3
"""Assertions for tests/check_resume_reference_change.sh.

Reads two Nextflow trace files -- a resume with nothing changed, and a resume after the
clb annotation's content changed -- and checks the F15 contract in both directions:

  * nothing changed        every task CACHED. Without this control, "the task re-ran"
                           below would be satisfied by a resume that caches nothing at all.
  * annotation changed     every process whose script carries the pks_annotation digest
                           re-ran, and the upstream processes that do not carry it are
                           still CACHED. A run that re-executed everything would also
                           satisfy "it re-ran", and would mean the digest is doing nothing.

Which processes carry the digest is read from the Modules/*.nf sources at run time, not
hard-coded here: adding the digest to another process must not require editing this file,
and removing it from one must make this test fail rather than quietly pass.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

# Processes that run before any annotation-dependent stage. They must stay cached, or the
# assertion below proves nothing about the digest -- only that resume is broken.
UPSTREAM_PROCESSES = ("filterReads", "mapReads")


def read_trace(path):
    """{process name: set(status)} from a Nextflow trace file."""
    statuses = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = (row.get("name") or "").split(" (")[0].strip()
            if not name:
                continue
            statuses.setdefault(name, set()).add((row.get("status") or "").strip())
    return statuses


def include_aliases(repo):
    """{declared process name: {names it is known by in the workflow}}.

    main.nf imports several processes under an alias -- `include { pksProfiler_align as
    pksProfilerAlign }` -- and the trace records the ALIAS. Matching the declared name
    against the trace would silently find nothing and report "this run executed no
    process carrying the digest", which looks like a broken test rather than the pass it
    actually is.
    """
    text = (Path(repo) / "main.nf").read_text()
    aliases = {}
    for block in re.finditer(r"include\s*\{([^}]*)\}\s*from", text):
        for item in block.group(1).split(";"):
            item = item.strip()
            if not item:
                continue
            parts = [p.strip() for p in item.split(" as ")]
            name = parts[0]
            aliases.setdefault(name, set()).add(parts[1] if len(parts) > 1 else name)
    return aliases


def processes_carrying(repo, digest):
    """Every process whose script block interpolates params.dep_digest?.<digest>, under
    every name the workflow knows it by."""
    declared = set()
    needle = f"dep_digest?.{digest}"
    for module in sorted((Path(repo) / "Modules").glob("*.nf")):
        text = module.read_text()
        for match in re.finditer(r"^process\s+(\w+)\s*\{", text, re.M):
            body = text[match.end():]
            end = body.find("\nprocess ")
            body = body[:end] if end > 0 else body
            if needle in body:
                declared.add(match.group(1))
    aliases = include_aliases(repo)
    carriers = set()
    for name in declared:
        carriers |= aliases.get(name, {name})
    return carriers


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--unchanged-trace", required=True)
    parser.add_argument("--changed-trace", required=True)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--digest", default="pks_annotation")
    args = parser.parse_args(argv)

    failures = []

    unchanged = read_trace(args.unchanged_trace)
    if not unchanged:
        failures.append(f"{args.unchanged_trace}: no tasks at all")
    not_cached = sorted(name for name, statuses in unchanged.items()
                        if statuses != {"CACHED"})
    if not_cached:
        failures.append(
            "resume with nothing changed re-ran: " + ", ".join(not_cached)
            + " -- every task should have been CACHED, so nothing below would be "
              "evidence about the dependency digest")

    carriers = processes_carrying(args.repo, args.digest)
    if not carriers:
        failures.append(f"no process in Modules/ carries the {args.digest} digest: "
                        f"the F15 mechanism this test is about is not present")

    changed = read_trace(args.changed_trace)
    if not changed:
        failures.append(f"{args.changed_trace}: no tasks at all")

    ran_again = sorted(name for name, statuses in changed.items() if statuses != {"CACHED"})
    still_cached = sorted(name for name, statuses in changed.items() if statuses == {"CACHED"})

    # Only carriers that this particular run actually executes can be checked.
    expected = sorted(carriers & set(changed))
    if not expected:
        failures.append(
            "this run executed no process carrying the digest "
            f"({', '.join(sorted(carriers))}); the check cannot conclude anything")
    for process in expected:
        if process not in ran_again:
            failures.append(
                f"{process} carries the {args.digest} digest but stayed CACHED after the "
                f"annotation's content changed -- resume reused work computed against the "
                f"old reference, which is exactly F15")

    for process in UPSTREAM_PROCESSES:
        if process in changed and process not in still_cached:
            failures.append(
                f"{process} re-ran although it does not depend on the annotation: the "
                f"digest is invalidating more than it should, so 're-ran' says nothing "
                f"about the dependency")

    print(f"processes carrying {args.digest}: {', '.join(sorted(carriers)) or 'none'}")
    print(f"re-ran after the reference change: {', '.join(ran_again) or 'none'}")
    print(f"still cached after the reference change: {', '.join(still_cached) or 'none'}")

    if failures:
        print("\nFAILURES:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print("\nOK: a reference content change invalidated exactly the tasks that depend on it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
