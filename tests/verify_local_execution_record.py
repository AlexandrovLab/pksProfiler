#!/usr/bin/env python3
"""Check the committed record of the execution checks a GitHub runner cannot run.

rc2 review item 2.2:

    Run real-execution tests on the actual release tag. The release CI should execute the
    end-to-end and mode-specific checks on the exact v1.0.0 tag, not only unit tests /
    graph construction.

Most of that is done on the runner: the tag build now runs check_e2e_execution.sh,
check_malformed_input_execution.sh, check_mode_execution.sh, check_v8_lanes_execution.sh
and check_resume_reference_change.sh on the exact tagged commit.

One check cannot go there. tests/check_metagenome_mag_execution.sh needs the real GTDB-Tk
r232 reference data (~100 GB), CheckM2 and geNomad; there is no small stand-in, and an
execution test whose entire subject is those databases being read by the real tools is
pointless against one. Fetching ~100 GB on every tag build is not a reasonable cost
either. The same applies to a zero-skip unit run: the GitHub runner has no R, so
tests/assert_required_coverage.py can only reach its budget there, never zero.

So, for those, the release gate requires a RECORDED result instead -- and verifies it
rather than trusting it:

  * the record carries a `source_digest`: a content digest of everything that can change
    what a run does (main.nf, nextflow.config, lib/, Modules/, scripts/, conf/,
    conda_envs/, ref/ and tests/, excluding the record and its logs). That digest is
    recomputed from the tree being released and must match. This is the assertion that
    matters, and it is deliberately not the commit sha: committing the record would
    change the commit, so a commit-based check could never be satisfied without an
    amend dance, while the source digest is unchanged by committing the record itself.
    It is also strictly stronger -- it is false for a dirty tree, which a sha is not;
  * the commit is recorded too, and reported, as human context for which commit the
    person was on;
  * every required check must be present with status `pass`;
  * each row carries the number of assertions/tests that ran, and a count of zero is
    refused -- a check that ran and asserted nothing is not evidence;
  * the record must be signed off by a person, by name. Producing it is a deliberate act
    on TSCC (tests/record_local_execution_checks.sh), not something CI can fake.

Usage
    python tests/verify_local_execution_record.py [--record FILE] [--source-digest D]
    python tests/verify_local_execution_record.py --print-source-digest

Exits non-zero, naming every problem, if the record does not stand up.
"""
import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECORD = ROOT / "tests/LOCAL_EXECUTION_CHECKS.tsv"

# Everything that can change what a run does. Deliberately NOT indices/ (binary Bowtie2
# index, unchanged since v0.0.1 and covered by params.dep_digest at run time) and not
# docs/ or the record's own logs.
SOURCE_PATHS = ("main.nf", "nextflow.config", "lib", "Modules", "scripts", "conf",
                "conda_envs", "ref", "tests")
SOURCE_EXCLUDE = {"tests/LOCAL_EXECUTION_CHECKS.tsv"}


def source_digest(root=ROOT):
    """sha256 over `<relative path>\\t<sha256 of content>` for every source file, sorted."""
    entries = []
    for name in SOURCE_PATHS:
        path = root / name
        if path.is_file():
            candidates = [path]
        elif path.is_dir():
            candidates = sorted(p for p in path.rglob("*") if p.is_file())
        else:
            continue
        for candidate in candidates:
            relative = candidate.relative_to(root).as_posix()
            if relative in SOURCE_EXCLUDE or "/__pycache__/" in f"/{relative}":
                continue
            digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
            entries.append(f"{relative}\t{digest}")
    return hashlib.sha256("\n".join(sorted(entries)).encode()).hexdigest()

# What the record must cover, and why it cannot be on the runner instead.
REQUIRED_CHECKS = {
    "tests/check_metagenome_mag_execution.sh":
        "needs GTDB-Tk r232 (~100 GB), CheckM2 and geNomad; no stand-in database can "
        "simulate the tool behaviour this check exists to exercise",
    "tests/assert_required_coverage.py":
        "a zero-skip unit run, which needs every real tool including R; the GitHub "
        "runner has no R and can only reach the skip budget, never zero",
}

HEADER_KEYS = ("commit", "source_digest", "recorded_on", "recorded_by", "host", "signed_off")
COLUMNS = ["check", "status", "tests", "duration_s", "log"]
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")


def parse(path):
    """({header key: value}, [{column: value}]) from the record."""
    header, rows, seen_columns = {}, [], None
    for number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            body = stripped.lstrip("#").strip()
            if ":" in body:
                key, _, value = body.partition(":")
                key = key.strip().lower().replace(" ", "_")
                if key in HEADER_KEYS:
                    header[key] = value.strip()
            continue
        fields = line.split("\t")
        if seen_columns is None:
            seen_columns = [f.strip() for f in fields]
            if seen_columns != COLUMNS:
                raise SystemExit(f"{path}:{number}: columns are {seen_columns}, "
                                 f"expected {COLUMNS}")
            continue
        if len(fields) != len(COLUMNS):
            raise SystemExit(f"{path}:{number}: {len(fields)} fields, expected {len(COLUMNS)}")
        rows.append(dict(zip(COLUMNS, (f.strip() for f in fields))))
    if seen_columns is None:
        raise SystemExit(f"{path}: no column header row")
    return header, rows


def head_commit():
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def verify(header, rows, digest):
    problems = []

    for key in HEADER_KEYS:
        if not header.get(key):
            problems.append(f"the record has no `{key}:` header")

    recorded_commit = header.get("commit", "")
    if recorded_commit and not COMMIT_RE.match(recorded_commit):
        problems.append(f"commit {recorded_commit!r} is not a 40-character sha")

    recorded_digest = header.get("source_digest", "")
    if recorded_digest and not re.match(r"^[0-9a-f]{64}$", recorded_digest):
        problems.append(f"source_digest {recorded_digest!r} is not a sha256")
    elif recorded_digest and recorded_digest != digest:
        problems.append(
            f"the record was produced from source {recorded_digest[:12]}..., but the "
            f"source being released digests to {digest[:12]}.... Something under "
            f"main.nf / lib / Modules / scripts / conf / conda_envs / ref / tests has "
            f"changed since these checks ran. Re-run "
            f"tests/record_local_execution_checks.sh on the tree being tagged.")

    signed = header.get("signed_off", "")
    if signed and len(signed) < 8:
        problems.append(f"signed_off is not a sign-off: {signed!r}")

    by_check = {}
    for row in rows:
        if row["check"] in by_check:
            problems.append(f"{row['check']} appears twice")
        by_check[row["check"]] = row

    for check, why in sorted(REQUIRED_CHECKS.items()):
        row = by_check.get(check)
        if row is None:
            problems.append(f"required check missing from the record: {check} ({why})")
            continue
        if row["status"] != "pass":
            problems.append(f"{check} is recorded as {row['status']!r}, not 'pass'")
        if not row["tests"].isdigit() or int(row["tests"]) <= 0:
            problems.append(f"{check} records {row['tests']!r} tests; a check that "
                            f"asserted nothing is not evidence")
        if not row["log"]:
            problems.append(f"{check} names no log")

    for row in rows:
        if row["status"] not in ("pass", "fail", "skipped"):
            problems.append(f"{row['check']}: unknown status {row['status']!r}")
        if row["status"] == "fail":
            problems.append(f"{row['check']} is recorded as FAILED; this is not releasable")

    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--record", default=str(DEFAULT_RECORD))
    parser.add_argument("--source-digest", default=None,
                        help="the digest to check against (default: computed from the tree)")
    parser.add_argument("--print-source-digest", action="store_true",
                        help="print the current tree's source digest and exit")
    args = parser.parse_args(argv)

    if args.print_source_digest:
        print(source_digest())
        return 0

    path = Path(args.record)
    if not path.is_file():
        print(f"no local execution-check record at {path}.\n"
              f"Produce one on a machine with the databases:\n"
              f"    bash tests/record_local_execution_checks.sh", file=sys.stderr)
        return 1

    header, rows = parse(path)
    digest = args.source_digest if args.source_digest is not None else source_digest()
    problems = verify(header, rows, digest)

    print(f"record         {path}")
    print(f"recorded from  source {header.get('source_digest', '(none)')}")
    print(f"releasing      source {digest}")
    here = head_commit()
    print(f"commit noted   {header.get('commit', '(none)')} "
          f"(context; the digest above is the check)")
    print(f"commit here    {here or '(not a git work tree)'}")
    print(f"recorded on    {header.get('recorded_on', '(none)')} by "
          f"{header.get('recorded_by', '(none)')} on {header.get('host', '(none)')}")
    print(f"signed off     {header.get('signed_off', '(none)')}")
    print()
    print(f"{'check':46s}{'status':10s}{'tests':>8s}{'seconds':>10s}")
    for row in rows:
        print(f"{row['check']:46s}{row['status']:10s}{row['tests']:>8s}{row['duration_s']:>10s}")

    if problems:
        print("\nLOCAL EXECUTION RECORD REJECTED:")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    print("\nLocal execution-check record accepted.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
