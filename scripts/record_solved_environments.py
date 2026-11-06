#!/usr/bin/env python3
"""Capture the SOLVED environments: the full transitive package list of each conda env the
release checks actually built, as reconstructable lock files.

rc2 review item 2.6 asked for pinned environments and, ideally, the solved environment or a
container digest. The ymls under conda_envs/ pin what we chose; this records what the solver
actually installed -- every package, build and channel URL with its md5 -- so an environment
can be recreated byte-for-byte with `conda create --file <lock>` even after channels move on.
No container path exists in this pipeline, and a container digest for an execution path the
validation never ran would be a claim, not a record; the lock files cover the path that ran.

A cache env carries no name, only a content hash, so each is mapped to its source yml by
requiring every pinned `name=version` spec of the yml to be installed in the env, taking the
yml with the LARGEST fully-satisfied pin set (pks_align_env pins a superset of minimap2_env;
the maximal rule keeps subset ymls from claiming the wrong env). Ambiguities and unmapped
envs are listed in INDEX.tsv rather than guessed.

    python3 scripts/record_solved_environments.py --cache <CONDA_CACHE_DIR> [--conda <bin>]

Writes conda_envs/solved/<yml stem>.lock.txt and conda_envs/solved/INDEX.tsv.
Ymls whose env is not in the cache are listed as not_captured; that is a fact about which
checks built which envs, not an error.
"""
import argparse
import hashlib
import json
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def yml_pins(path):
    """{package: version} for every `- name=version` spec in a conda yml."""
    pins = {}
    for line in path.read_text().splitlines():
        m = re.match(r"\s*-\s*([A-Za-z0-9_.-]+)={1,2}([A-Za-z0-9_.]+)\s*(#.*)?$", line)
        if m:
            pins[m.group(1).lower()] = m.group(2)
    return pins


def installed(env_dir):
    """{package: version} from conda-meta/*.json, the ground truth of what is in an env."""
    out = {}
    for meta in (env_dir / "conda-meta").glob("*.json"):
        try:
            rec = json.loads(meta.read_text())
            out[rec["name"].lower()] = rec["version"]
        except (ValueError, KeyError):
            pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--conda", default="conda")
    a = ap.parse_args()

    cache = Path(a.cache)
    envs = [d for d in sorted(cache.iterdir())
            if d.is_dir() and (d / "conda-meta" / "history").exists()]
    ymls = {p.stem: yml_pins(p) for p in sorted((ROOT / "conda_envs").glob("*.yml"))}

    out_dir = ROOT / "conda_envs" / "solved"
    out_dir.mkdir(exist_ok=True)
    rows, claimed = [], {}

    for env in envs:
        have = installed(env)
        full = {name: pins for name, pins in ymls.items()
                if pins and all(have.get(p) == v for p, v in pins.items())}
        if not full:
            rows.append(("-", env.name, "unmapped", "-", "no yml's pins are all present"))
            continue
        best_size = max(len(p) for p in full.values())
        best = [n for n, p in full.items() if len(p) == best_size]
        if len(best) > 1:
            rows.append(("-", env.name, "ambiguous", "-", "ties: " + ",".join(sorted(best))))
            continue
        name = best[0]
        if name in claimed:
            rows.append(("-", env.name, "duplicate", "-", f"{name} already mapped to {claimed[name]}"))
            continue
        claimed[name] = env.name

        lock = subprocess.run([a.conda, "list", "--explicit", "--md5", "-p", str(env)],
                              capture_output=True, text=True, check=True).stdout
        n_pkgs = sum(1 for l in lock.splitlines() if l.startswith("https://"))
        header = (f"# Solved environment for conda_envs/{name}.yml\n"
                  f"# Captured {date.today().isoformat()} from the frozen release-check cache\n"
                  f"# (cache env {env.name}; {n_pkgs} packages).\n"
                  f"# Recreate:  conda create --prefix <dir> --file <this file>\n")
        path = out_dir / f"{name}.lock.txt"
        path.write_text(header + lock)
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append((name, env.name, "captured", sha, f"{n_pkgs} packages"))
        print(f"  {name:28s} <- {env.name}  ({n_pkgs} packages)")

    for name in sorted(set(ymls) - set(claimed)):
        rows.append((name, "-", "not_captured", "-",
                     "its env is not in this cache (not built by the recorded checks)"))

    index = out_dir / "INDEX.tsv"
    with index.open("w") as fh:
        fh.write(f"# Written by scripts/record_solved_environments.py on {date.today().isoformat()}; "
                 f"cache: {cache}\n")
        fh.write("yml\tcache_env\tstatus\tlock_sha256\tnote\n")
        for r in sorted(rows):
            fh.write("\t".join(r) + "\n")
    captured = sum(1 for r in rows if r[2] == "captured")
    print(f"{captured} captured, {sum(1 for r in rows if r[2]=='not_captured')} not in cache "
          f"-> {index}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
