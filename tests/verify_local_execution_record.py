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
  * every required check must be present with status `pass` -- or, for a check other
    than the unit gate itself, status `skipped` together with a `skip_waiver:` header
    line naming the check and giving a substantive, signed-off reason. The waiver is
    part of the record, so it is as visible and as attributable as a pass;
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

# The unit gate can never be waived: it runs anywhere Python runs, so there is no
# environment that justifies skipping it, and it is the check that counts the skips.
NEVER_SKIPPABLE = {"tests/assert_required_coverage.py"}
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
                elif key == "skip_waiver":
                    # "# skip_waiver: <check>: [since=<sha>] <reason>", once per check.
                    #
                    # since= is what stops a rationale outliving the tree it was reasoned
                    # about. The v1.0.0 record carried a waiver saying nothing the MAG
                    # check exercises had changed since 2026-10-22; by then the tree held
                    # a later fix to Modules/pks_mag.nf, which is the module that check
                    # exists to exercise, and Modules/ is inside the source digest. The
                    # sentence was false and nothing could tell. With since=, the verifier
                    # recomputes what changed and expires the waiver itself.
                    check, _, reason = value.strip().partition(":")
                    reason = reason.strip()
                    since = ""
                    m = re.match(r"since=(\S+)\s+(.*)", reason)
                    if m:
                        since, reason = m.group(1), m.group(2).strip()
                    header.setdefault("skip_waivers", {})[check.strip()] = {
                        "reason": reason, "since": since}
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


def digest_scope_changes_since(commit, root=ROOT):
    """Files inside SOURCE_PATHS that differ between `commit` and the working tree.

    None if the commit cannot be resolved. The record itself and its logs are excluded
    for the same reason the digest excludes them: committing the record must not expire
    the waiver that record carries.
    """
    try:
        done = subprocess.run(["git", "-C", str(root), "diff", "--name-only", commit, "--"]
                              + list(SOURCE_PATHS),
                              capture_output=True, text=True)
    except OSError:
        return None
    if done.returncode != 0:
        return None
    changed = set()
    for name in done.stdout.split("\n"):
        name = name.strip()
        if not name or any(name == x or name.startswith(x.rstrip("/") + "/")
                           for x in SOURCE_EXCLUDE):
            continue
        changed.add(name)
    return changed


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

    waivers = header.get("skip_waivers", {})
    for check, why in sorted(REQUIRED_CHECKS.items()):
        row = by_check.get(check)
        if row is None:
            problems.append(f"required check missing from the record: {check} ({why})")
            continue
        if row["status"] == "skipped":
            waiver = waivers.get(check) or {}
            if check in NEVER_SKIPPABLE:
                problems.append(f"{check} may not be skipped: it is the unit gate itself")
            elif len(waiver.get("reason", "")) < 20:
                problems.append(
                    f"{check} is recorded as skipped without a skip_waiver header naming "
                    f"it and giving a substantive reason; a skip nobody wrote down is "
                    f"not a waiver")
            elif not waiver.get("since"):
                problems.append(
                    f"{check}'s waiver does not say which tree it was reasoned against. "
                    f"Add `since=<commit>` after the check name, naming the commit whose "
                    f"state the reason describes, so a later change to the code this "
                    f"check exercises expires the waiver instead of outliving it")
            else:
                changed = digest_scope_changes_since(waiver["since"])
                if changed is None:
                    problems.append(
                        f"{check}'s waiver names since={waiver['since']}, which is not a "
                        f"commit in this repository")
                elif changed:
                    shown = ", ".join(sorted(changed)[:5])
                    more = f" (+{len(changed) - 5} more)" if len(changed) > 5 else ""
                    problems.append(
                        f"{check}'s waiver was reasoned against {waiver['since']}, but "
                        f"{len(changed)} file(s) inside the source digest have changed "
                        f"since then: {shown}{more}. The reason may no longer be true, so "
                        f"it has expired: re-run the check, or write a waiver that "
                        f"describes THIS tree")
            continue  # a waived skip asserts nothing, so the tests/log demands below do not apply
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
    for check, waiver in sorted(header.get("skip_waivers", {}).items()):
        since = waiver.get("since") or "(no tree named)"
        print(f"WAIVED        {check}   reasoned against {since}")
        print(f"              {waiver.get('reason', '')}")
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
