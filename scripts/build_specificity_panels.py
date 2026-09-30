#!/usr/bin/env python3
"""Build the two pks specificity reference panels from public sequence, and screen every
member against the canonical clb island before it is allowed into a panel.

Ludmil, 25 Sep report, section 4, asks two different questions that were previously one
"specificity check", and says explicitly that they must not be conflated:

  V8a  LOCUS-DETECTION SPECIFICITY -- "does this sequence really belong to the pks/clb
       locus?"  Challenge the supporting reads against pks-NEGATIVE Enterobacterales,
       close homologous NRPS/PKS sequences, broad gut/oral backgrounds, GRCh38 and
       UniVec.  Panel: ref/specificity_panels/v8a_locus_decoy_panel.tsv.

  V8b  CARRIER ATTRIBUTION -- "which organism is the likely carrier?"  Compare the same
       reads against alternative pks-POSITIVE carriers and strains (pks+ E. coli beyond
       IHE3034, pks+ Klebsiella, Citrobacter, Enterobacter).  Panel:
       ref/specificity_panels/v8b_carrier_panel.tsv.

The two panels have opposite membership tests, and this script applies them rather than
trusting the accession lists:

    island coverage of a member  =  fraction of the canonical 50,768 bp IHE3034 island
                                    (params.pks_start_1based..pks_end_1based) covered by
                                    minimap2 -cx asm20 alignments of the island against
                                    that member's assembly

    V8a decoy   accepted only if coverage <  --max-decoy-island-coverage   (default 0.10)
    V8b carrier accepted only if coverage >= --min-carrier-island-coverage (default 0.50)

A rejected member is never silently dropped: it is written to <panel>_pks_screen.tsv with
its measured coverage and a status, and left out of the FASTA and the index. Both
directions matter. A pks-carrying "decoy" absorbs genuine pks reads and turns V8a into a
coin flip; a pks-negative "carrier" can never win a read and only inflates V8b.

Everything here is public reference sequence -- NCBI Datasets assemblies, NCBI UniVec,
MIBiG biosynthetic gene clusters, and (optionally) a local GRCh38 FASTA. No sample data
is read, produced or touched.

    # V8a: pks-negative decoys, NRPS/PKS homolog genomes, gut/oral backgrounds, UniVec
    build_specificity_panels.py --panel v8a \\
        --members ref/specificity_panels/v8a_locus_decoy_panel.tsv \\
        --island-fasta island.fna --outdir .../v8a --minimap2 $(command -v minimap2) \\
        --univec [--grch38 /path/GRCh38.fa] [--mibig-gbk mibig_gbk_4.0.tar.gz]

    # V8b: pks-positive carriers
    build_specificity_panels.py --panel v8b \\
        --members ref/specificity_panels/v8b_carrier_panel.tsv \\
        --island-fasta island.fna --outdir .../v8b --minimap2 $(command -v minimap2)

Outputs, per panel (all paths relative to --outdir):

  members/<key>.fna              one downloaded assembly per accepted member, headers
                                 retagged '>{key}__{original}' so every contig names its
                                 own member for the rest of the pipeline
  <panel>_panel.fna              the accepted members concatenated
  <panel>_panel.mmi              minimap2 -x sr index of that FASTA
  <panel>_contigs.txt            one contig name per line (V8a: the decoy-side contig
                                 list competitive_realign_score.py reads)
  <panel>_carrier_map.tsv        V8b only: contig -> carrier key/organism/accession
  <panel>_panel_members.tsv      key, accession, category, organism, rationale, bases,
                                 contigs, md5, status
  <panel>_pks_screen.tsv         key, island_coverage, aligned_bp, status -- every
                                 candidate, accepted and rejected alike
  README.md                      what was built, from what, when, and with which tools

Concatenated true+decoy index (V8a only, with --true-fasta): the two-index Arm B in
scripts/competitive_realign_score.py also needs one index holding the true reference and
the decoys together, plus the matching true-contig list. Passing --true-fasta writes
v8a_locus_concatenated.fna/.mmi, v8a_true_contigs.txt and v8a_decoy_contigs.txt as well.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"
UNIVEC_URL = "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"
MIBIG_DEFAULT_URL = "https://dl.secondarymetabolites.org/mibig/mibig_gbk_4.0.tar.gz"

# MIBiG clusters whose product is colibactin itself: the one thing that must never enter
# a decoy panel built to challenge colibactin reads.
COLIBACTIN_PATTERN = re.compile(r"colibactin|\bclb[A-S]\b", re.IGNORECASE)
# The biosynthetic classes V8a is about: NRPS and PKS megasynthase machinery whose
# domains clbB/clbC/clbH/clbI/clbJ/clbK/clbN/clbO share with unrelated clusters.
NRPS_PKS_PATTERN = re.compile(r"\bNRPS\b|nonribosomal|polyketide|\bPKS\b|\bT1PKS\b|transAT",
                              re.IGNORECASE)


# --------------------------------------------------------------------------- utilities
def sh(cmd, **kwargs):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kwargs)


def md5_of(path, chunk=1 << 20):
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url, dest, retries=4, sleep=3):
    """Download to dest, retrying: NCBI rate-limits, and a half-written file here would
    silently become a truncated panel member."""
    last = None
    for attempt in range(retries):
        try:
            tmp = Path(str(dest) + ".part")
            with urllib.request.urlopen(url, timeout=300) as response, tmp.open("wb") as out:
                shutil.copyfileobj(response, out)
            tmp.replace(dest)
            return dest
        except Exception as exc:  # noqa: BLE001 - reported, then retried or raised
            last = exc
            time.sleep(sleep * (attempt + 1))
    raise RuntimeError(f"could not download {url}: {last}")


def read_members(path):
    rows = []
    with Path(path).open() as handle:
        lines = [ln for ln in handle if not ln.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        if row.get("key"):
            rows.append(row)
    if not rows:
        raise SystemExit(f"no member rows in {path}")
    return rows


def read_fasta(path):
    """Yield (header, sequence) for a possibly large FASTA, without holding all of it."""
    header, seq = None, []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def write_tagged(src, dest, key):
    """Copy a FASTA, prefixing every contig name with '<key>__' so a contig always names
    the panel member it came from. Returns (n_contigs, n_bases, [contig names])."""
    contigs, bases = [], 0
    with Path(dest).open("w") as out:
        for header, seq in read_fasta(src):
            name = f"{key}__{header.split()[0]}"
            contigs.append(name)
            bases += len(seq)
            out.write(f">{name} {header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
    return len(contigs), bases, contigs


# ------------------------------------------------------------------------- NCBI genomes
def assembly_report(accession):
    url = f"{DATASETS_API}/{urllib.parse.quote(accession)}/dataset_report"
    with urllib.request.urlopen(url, timeout=120) as response:
        payload = json.load(response)
    reports = payload.get("reports") or []
    if not reports:
        raise RuntimeError(f"NCBI Datasets returned no report for {accession}")
    report = reports[0]
    return {
        "accession": report.get("accession"),
        "organism": report.get("organism", {}).get("organism_name", ""),
        "strain": (report.get("organism", {}).get("infraspecific_names") or {}).get("strain", ""),
        "assembly_name": report.get("assembly_info", {}).get("assembly_name", ""),
        "assembly_level": report.get("assembly_info", {}).get("assembly_level", ""),
        "assembly_status": report.get("assembly_info", {}).get("assembly_status", ""),
        "total_length": report.get("assembly_stats", {}).get("total_sequence_length", ""),
    }


def download_assembly(accession, dest_fna, cache_dir):
    """Fetch GENOME_FASTA for one accession through the NCBI Datasets v2 download
    endpoint and unpack the single .fna it contains."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    zip_path = cache_dir / f"{accession}.zip"
    if not zip_path.exists() or zip_path.stat().st_size == 0:
        url = (f"{DATASETS_API}/{urllib.parse.quote(accession)}/download"
               "?include_annotation_type=GENOME_FASTA&filename=genome.zip")
        fetch(url, zip_path)
    with zipfile.ZipFile(zip_path) as archive:
        names = [n for n in archive.namelist()
                 if n.endswith((".fna", ".fa", ".fasta")) and "/data/" in n]
        if not names:
            raise RuntimeError(f"{accession}: no FASTA inside the Datasets zip")
        names.sort()
        with archive.open(names[0]) as src, Path(dest_fna).open("wb") as out:
            shutil.copyfileobj(src, out)
    return dest_fna


# ------------------------------------------------------------------------ island screen
def island_coverage(minimap2, island_fasta, genome_fasta, preset, threads):
    """Fraction of the island covered by alignments of the island against this assembly.

    Query is the island, target is the candidate genome, so PAF columns 3/4 are island
    coordinates and the union of those intervals is what "how much of the locus does this
    organism carry" means. Overlapping alignments are merged, never summed.
    """
    result = subprocess.run(
        [minimap2, "-c", "-x", preset, "-t", str(threads), "--secondary=no",
         str(genome_fasta), str(island_fasta)],
        check=True, capture_output=True, text=True)
    island_len = sum(len(seq) for _h, seq in read_fasta(island_fasta))
    intervals = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        intervals.append((int(fields[2]), int(fields[3])))
    intervals.sort()
    covered, merged_end = 0, 0
    for start, end in intervals:
        if end <= merged_end:
            continue                      # wholly inside an interval already counted
        covered += end - max(start, merged_end)
        merged_end = end
    return (covered / island_len if island_len else 0.0), covered, island_len


# ------------------------------------------------------------------------------- MIBiG
def mibig_clusters(tarball, out_dir, max_clusters):
    """Extract non-colibactin NRPS/PKS biosynthetic gene clusters from a MIBiG GenBank
    release as one FASTA per cluster.

    MIBiG ships one GenBank record per BGC. A cluster is kept when its text names an NRPS
    or PKS biosynthetic class and does NOT name colibactin or a clb gene -- the colibactin
    cluster itself is the one thing that cannot be a decoy for colibactin reads.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    kept, skipped_colibactin, seen = [], 0, 0
    with tarfile.open(tarball, "r:gz") as archive:
        # The whole archive is scanned even once --mibig-max clusters are in hand: the
        # colibactin cluster could sit anywhere in it, and "0 colibactin clusters
        # excluded" must mean the archive holds none that matched, not that the scan
        # stopped before reaching one. Tar order is fixed, so the selection is
        # reproducible.
        for entry in archive:
            if not entry.isfile() or not entry.name.endswith((".gbk", ".gb")):
                continue
            seen += 1
            handle = archive.extractfile(entry)
            if handle is None:
                continue
            text = handle.read().decode("utf-8", errors="replace")
            bgc_id = Path(entry.name).stem.split(".")[0]
            if COLIBACTIN_PATTERN.search(text):
                skipped_colibactin += 1
                continue
            if not NRPS_PKS_PATTERN.search(text) or len(kept) >= max_clusters:
                continue
            sequence = _genbank_sequence(text)
            if len(sequence) < 5000:
                continue
            path = out_dir / f"{bgc_id}.fna"
            with path.open("w") as out:
                out.write(f">{bgc_id} MIBiG non-colibactin NRPS/PKS cluster\n")
                for i in range(0, len(sequence), 60):
                    out.write(sequence[i:i + 60] + "\n")
            kept.append((bgc_id, path, len(sequence)))
    return kept, seen, skipped_colibactin


def _genbank_sequence(text):
    if "ORIGIN" not in text:
        return ""
    body = text.split("ORIGIN", 1)[1]
    body = body.split("//", 1)[0]
    return re.sub(r"[^acgtnACGTN]", "", body).upper()


# -------------------------------------------------------------------------------- build
def build(args):
    out = Path(args.outdir)
    (out / "members").mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache_dir) if args.cache_dir else out / "_downloads"

    members = read_members(args.members)
    accepted, screen_rows, member_rows, all_contigs = [], [], [], []

    for row in members:
        key, accession = row["key"], row["accession"]
        raw = cache / f"{key}.raw.fna"
        if not raw.exists() or raw.stat().st_size == 0:
            raw.parent.mkdir(parents=True, exist_ok=True)
            print(f"[download] {key} {accession}", flush=True)
            download_assembly(accession, raw, cache)
        report = assembly_report(accession)
        expected = (row.get("expected_organism") or "").strip()
        organism = report["organism"]
        # Phase-1 convention, kept: identity is confirmed against the assembly's own
        # report, never assumed from the accession string.
        genus_expected = expected.split()[0] if expected else ""
        identity_ok = bool(genus_expected) and organism.startswith(genus_expected)

        coverage, aligned_bp, island_len = island_coverage(
            args.minimap2, args.island_fasta, raw, args.screen_preset, args.threads)

        if args.panel == "v8a":
            ok = coverage < args.max_decoy_island_coverage
            status = "accepted_pks_negative_decoy" if ok else "rejected_pks_positive"
        else:
            ok = coverage >= args.min_carrier_island_coverage
            status = "accepted_pks_positive_carrier" if ok else "rejected_not_pks_positive"
        if not identity_ok:
            ok, status = False, "rejected_organism_mismatch"

        screen_rows.append({
            "key": key, "accession": accession, "organism": organism,
            "expected_organism": expected, "identity_confirmed": identity_ok,
            "island_bp_covered": aligned_bp, "island_bp_total": island_len,
            "island_coverage": f"{coverage:.6f}", "screen_preset": args.screen_preset,
            "status": status,
        })
        print(f"[screen]   {key:<40} coverage={coverage:.4f} {status}", flush=True)
        if not ok:
            continue

        tagged = out / "members" / f"{key}.fna"
        n_contigs, n_bases, contigs = write_tagged(raw, tagged, key)
        all_contigs += contigs
        accepted.append(tagged)
        member_rows.append({
            "key": key, "accession": accession,
            "category": row.get("category") or row.get("genus_group", ""),
            "organism": organism, "strain": report["strain"],
            "assembly_name": report["assembly_name"],
            "assembly_level": report["assembly_level"],
            "contigs": n_contigs, "bases": n_bases, "md5": md5_of(tagged),
            "island_coverage": f"{coverage:.6f}", "status": status,
            "rationale": row.get("rationale", ""),
        })

    # --- V8a-only extra members: UniVec, GRCh38, MIBiG clusters -----------------------
    if args.panel == "v8a":
        if args.univec:
            raw = cache / "univec.fna"
            if not raw.exists() or raw.stat().st_size == 0:
                print("[download] UniVec", flush=True)
                fetch(UNIVEC_URL, raw)
            tagged = out / "members" / "univec.fna"
            n_contigs, n_bases, contigs = write_tagged(raw, tagged, "univec")
            all_contigs += contigs
            accepted.append(tagged)
            member_rows.append({
                "key": "univec", "accession": "NCBI UniVec", "category": "vector",
                "organism": "vector/adapter/linker sequences", "strain": "",
                "assembly_name": "UniVec", "assembly_level": "",
                "contigs": n_contigs, "bases": n_bases, "md5": md5_of(tagged),
                "island_coverage": "", "status": "accepted_pks_negative_decoy",
                "rationale": "Ludmil V8a names UniVec explicitly: a pks read must not be "
                             "explainable as vector, adapter or linker sequence.",
            })
        if args.grch38:
            tagged = out / "members" / "grch38.fna"
            print(f"[tag]      GRCh38 from {args.grch38} (this is the slow one)", flush=True)
            n_contigs, n_bases, contigs = write_tagged(args.grch38, tagged, "grch38")
            all_contigs += contigs
            accepted.append(tagged)
            member_rows.append({
                "key": "grch38", "accession": str(args.grch38), "category": "host",
                "organism": "Homo sapiens GRCh38", "strain": "",
                "assembly_name": Path(args.grch38).name, "assembly_level": "",
                "contigs": n_contigs, "bases": n_bases, "md5": md5_of(tagged),
                "island_coverage": "", "status": "accepted_pks_negative_decoy",
                "rationale": "Ludmil V8a names GRCh38 explicitly: the host genome the "
                             "cohort was aligned against is the largest single "
                             "opportunity for a spurious pks read.",
            })
        if args.mibig_gbk:
            # Path() collapses "https://" to "https:/", so the URL test happens on the
            # raw string, before it is ever turned into a path.
            if args.mibig_gbk.startswith(("http://", "https://")):
                tarball = cache / "mibig_gbk.tar.gz"
                if not tarball.exists() or tarball.stat().st_size == 0:
                    print(f"[download] MIBiG {args.mibig_gbk}", flush=True)
                    fetch(args.mibig_gbk, tarball)
            else:
                tarball = Path(args.mibig_gbk)
            clusters, seen, skipped = mibig_clusters(
                tarball, cache / "mibig_clusters", args.mibig_max)
            print(f"[mibig]    {len(clusters)} NRPS/PKS clusters kept of {seen} "
                  f"({skipped} colibactin clusters excluded)", flush=True)
            combined = out / "members" / "mibig_nrps_pks.fna"
            contigs, bases = [], 0
            with combined.open("w") as handle:
                for bgc_id, path, length in clusters:
                    for header, seq in read_fasta(path):
                        name = f"mibig__{header.split()[0]}"
                        contigs.append(name)
                        bases += len(seq)
                        handle.write(f">{name} {header}\n")
                        for i in range(0, len(seq), 60):
                            handle.write(seq[i:i + 60] + "\n")
            if contigs:
                all_contigs += contigs
                accepted.append(combined)
                member_rows.append({
                    "key": "mibig_nrps_pks", "accession": str(args.mibig_gbk),
                    "category": "nrps_pks_homolog_cluster",
                    "organism": f"MIBiG non-colibactin NRPS/PKS clusters (n={len(clusters)})",
                    "strain": "", "assembly_name": "MIBiG", "assembly_level": "",
                    "contigs": len(contigs), "bases": bases, "md5": md5_of(combined),
                    "island_coverage": "", "status": "accepted_pks_negative_decoy",
                    "rationale": "Ludmil V8a names close homologous NRPS/PKS sequences: "
                                 "the megasynthase domains clbB/clbC/clbH/clbI/clbJ/clbK/"
                                 "clbN/clbO share with unrelated clusters, with the "
                                 "colibactin cluster itself excluded by construction.",
                })

    if not accepted:
        raise SystemExit("every candidate failed the island screen; no panel was built")

    # --- concatenate, index, write the tables ----------------------------------------
    prefix = args.panel
    panel_fna = out / f"{prefix}_panel.fna"
    with panel_fna.open("wb") as out_handle:
        for path in accepted:
            with Path(path).open("rb") as src:
                shutil.copyfileobj(src, out_handle)

    contigs_txt = out / f"{prefix}_contigs.txt"
    contigs_txt.write_text("\n".join(all_contigs) + "\n")

    if not args.no_index:
        print(f"[index]    minimap2 -x sr -d {panel_fna.name}", flush=True)
        sh([args.minimap2, "-x", "sr", "-t", str(args.threads),
            "-d", str(out / f"{prefix}_panel.mmi"), str(panel_fna)])

    if args.panel == "v8b":
        with (out / f"{prefix}_carrier_map.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["contig", "carrier_key", "carrier_organism", "carrier_accession"])
            by_key = {row["key"]: row for row in member_rows}
            for contig in all_contigs:
                key = contig.split("__", 1)[0]
                row = by_key.get(key, {})
                writer.writerow([contig, key, row.get("organism", ""), row.get("accession", "")])

    if args.panel == "v8a" and args.true_fasta:
        # Arm B's concatenated architecture: true reference + decoys in ONE index, plus
        # the two contig lists that say which side of it each contig belongs to.
        true_tagged = out / "members" / "true_reference.fna"
        n_true, _bases, true_contigs = write_tagged(args.true_fasta, true_tagged, "true")
        concatenated = out / f"{prefix}_locus_concatenated.fna"
        with concatenated.open("wb") as out_handle:
            with true_tagged.open("rb") as src:
                shutil.copyfileobj(src, out_handle)
            with panel_fna.open("rb") as src:
                shutil.copyfileobj(src, out_handle)
        (out / f"{prefix}_true_contigs.txt").write_text("\n".join(true_contigs) + "\n")
        (out / f"{prefix}_decoy_contigs.txt").write_text("\n".join(all_contigs) + "\n")
        if not args.no_index:
            print(f"[index]    minimap2 -x sr -d {concatenated.name}", flush=True)
            sh([args.minimap2, "-x", "sr", "-t", str(args.threads),
                "-d", str(out / f"{prefix}_locus_concatenated.mmi"), str(concatenated)])
            sh([args.minimap2, "-x", "sr", "-t", str(args.threads),
                "-d", str(out / f"{prefix}_true_index.mmi"), str(true_tagged)])
        print(f"[true]     {n_true} true-reference contig(s)", flush=True)

    member_fields = ["key", "accession", "category", "organism", "strain", "assembly_name",
                     "assembly_level", "contigs", "bases", "md5", "island_coverage",
                     "status", "rationale"]
    with (out / f"{prefix}_panel_members.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, member_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(member_rows)

    screen_fields = ["key", "accession", "organism", "expected_organism",
                     "identity_confirmed", "island_bp_covered", "island_bp_total",
                     "island_coverage", "screen_preset", "status"]
    with (out / f"{prefix}_pks_screen.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, screen_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(screen_rows)

    write_readme(args, out, member_rows, screen_rows, panel_fna)
    rejected = [r for r in screen_rows if r["status"].startswith("rejected")]
    print(f"\n{prefix}: {len(member_rows)} members accepted, {len(rejected)} rejected, "
          f"{len(all_contigs)} contigs, {panel_fna.stat().st_size / 1e6:.1f} MB FASTA")
    print(f"panel: {panel_fna}")
    for row in rejected:
        print(f"  rejected {row['key']}: {row['status']} (island coverage {row['island_coverage']})")
    return 0


def write_readme(args, out, member_rows, screen_rows, panel_fna):
    minimap2_version = subprocess.run([args.minimap2, "--version"], capture_output=True,
                                       text=True).stdout.strip()
    title = ("V8a locus-detection specificity decoy panel"
             if args.panel == "v8a" else "V8b carrier-attribution panel")
    question = ("Does a pks-supporting read really belong to the pks/clb locus, or is it "
                "equally explainable by a pks-negative relative, a homologous NRPS/PKS "
                "cluster, a gut/oral background, the host genome or a vector?"
                if args.panel == "v8a" else
                "Given that a read does belong to the pks/clb locus, WHICH pks-positive "
                "organism is the likely carrier? This is not a locus-specificity test and "
                "must not be reported as one.")
    lines = [
        f"# {title}",
        "",
        f"Built {time.strftime('%Y-%m-%d %H:%M')} by `scripts/build_specificity_panels.py` "
        f"from `{Path(args.members).name}`.",
        "",
        f"**Question this panel answers.** {question}",
        "",
        "Every file here is public reference sequence (NCBI Datasets assemblies, NCBI "
        "UniVec, MIBiG clusters, and any local GRCh38 passed in). No sample-derived data.",
        "",
        "## Membership test",
        "",
        "Members are not trusted from the accession list: each assembly is screened by "
        f"aligning the canonical IHE3034 clb island against it (`minimap2 -cx "
        f"{args.screen_preset}`) and measuring the merged fraction of the island covered.",
        "",
        (f"* V8a decoy accepted when island coverage < {args.max_decoy_island_coverage}"
         if args.panel == "v8a" else
         f"* V8b carrier accepted when island coverage >= {args.min_carrier_island_coverage}"),
        "* organism identity confirmed against the assembly's own NCBI report",
        f"* full result for every candidate, accepted and rejected: `{args.panel}_pks_screen.tsv`",
        "",
        "## Accepted members",
        "",
        "| key | accession | category | organism | contigs | bases | island coverage |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in member_rows:
        lines.append(f"| {row['key']} | {row['accession']} | {row['category']} | "
                     f"{row['organism']} | {row['contigs']} | {row['bases']} | "
                     f"{row['island_coverage'] or 'n/a'} |")
    rejected = [r for r in screen_rows if r["status"].startswith("rejected")]
    lines += ["", "## Rejected candidates", ""]
    if rejected:
        lines += ["| key | accession | island coverage | why |", "|---|---|---|---|"]
        for row in rejected:
            lines.append(f"| {row['key']} | {row['accession']} | {row['island_coverage']} | "
                         f"{row['status']} |")
    else:
        lines.append("None: every candidate passed the screen.")
    lines += [
        "",
        "## Files",
        "",
        f"* `{args.panel}_panel.fna` / `.mmi` -- the panel and its minimap2 `-x sr` index "
        f"({panel_fna.stat().st_size / 1e6:.1f} MB FASTA)",
        f"* `{args.panel}_contigs.txt` -- one contig per line",
        (f"* `{args.panel}_carrier_map.tsv` -- contig to carrier, read by "
         "`scripts/carrier_attribution.py`" if args.panel == "v8b" else
         f"* `{args.panel}_locus_concatenated.fna`/`.mmi`, `{args.panel}_true_contigs.txt`, "
         f"`{args.panel}_decoy_contigs.txt` -- the concatenated Arm B architecture "
         "(written only when --true-fasta was given)"),
        f"* `{args.panel}_panel_members.tsv`, `{args.panel}_pks_screen.tsv` -- provenance",
        "",
        "## Tools",
        "",
        f"* minimap2 {minimap2_version}",
        f"* screen preset `{args.screen_preset}`, index preset `sr`",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines))


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel", required=True, choices=["v8a", "v8b"])
    p.add_argument("--members", required=True, help="panel TSV from ref/specificity_panels/")
    p.add_argument("--outdir", required=True)
    p.add_argument("--island-fasta", required=True,
                   help="the canonical clb island FASTA the screen aligns against every "
                        "candidate (samtools faidx <pks_reference_fasta> "
                        "NC_017628.1:2193827-2244594)")
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--cache-dir", default=None,
                   help="where raw downloads are kept (default <outdir>/_downloads); "
                        "reused across runs so a rebuild does not re-download")
    p.add_argument("--screen-preset", default="asm20",
                   help="minimap2 preset for the island screen: asm20 tolerates ~5%% "
                        "divergence, which is the right order for a mobile island across "
                        "genera (default: asm20)")
    p.add_argument("--max-decoy-island-coverage", type=float, default=0.10,
                   help="V8a: a decoy must cover LESS than this fraction of the island")
    p.add_argument("--min-carrier-island-coverage", type=float, default=0.50,
                   help="V8b: a carrier must cover AT LEAST this fraction of the island")
    p.add_argument("--univec", action="store_true", help="V8a: add NCBI UniVec")
    p.add_argument("--grch38", default=None,
                   help="V8a: add a local GRCh38 FASTA. Large -- roughly 3 GB of sequence "
                        "and a ~7 GB index -- so it is opt-in, not a default.")
    p.add_argument("--mibig-gbk", default=None,
                   help=f"V8a: MIBiG GenBank release tarball, path or URL "
                        f"(e.g. {MIBIG_DEFAULT_URL}). Non-colibactin NRPS/PKS clusters "
                        f"are extracted from it.")
    p.add_argument("--mibig-max", type=int, default=400,
                   help="cap on MIBiG clusters added (default: 400)")
    p.add_argument("--true-fasta", default=None,
                   help="V8a: the true pks+ reference FASTA; also writes the concatenated "
                        "true+decoy index and the two contig lists Arm B needs")
    p.add_argument("--no-index", action="store_true",
                   help="write the FASTAs and tables but skip minimap2 indexing")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.panel == "v8b" and (args.univec or args.grch38 or args.mibig_gbk):
        raise SystemExit("--univec/--grch38/--mibig-gbk belong to the V8a locus panel; "
                          "the V8b carrier panel holds pks-positive carriers only")
    return build(args)


if __name__ == "__main__":
    sys.exit(main())
