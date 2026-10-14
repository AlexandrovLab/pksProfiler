"""T7 (Ludmil, 28 Sep): the run records which database build it read.

"Freeze and record versions/checksums of KrakenUniq/Bracken, GTDB, CheckM2, geNomad,
MLST and any other DB."

Before this, RUN_REPORT's DATABASES block printed the PATH each database was read from.
A path is not a version: unpack GTDB r233 over .../gtdbtk_r232 and every run before and
after reports the same line, so a results tree cannot say which reference it used.

Three pieces, tested here without Nextflow:

  conf/databases.tsv                      the frozen manifest: release, size, file count,
                                          fingerprint, and a real content sha256 of one
                                          small key file per database
  scripts/record_database_versions.sh     recomputes the fingerprint at run time and
                                          compares it with the manifest
  finalize_provenance.render_database_versions  puts it in RUN_REPORT

The cases that matter: a database that matches the frozen manifest, one that has been
rebuilt (the whole point), one the manifest has never seen, an unset parameter, and a
path that is not there at all.
"""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/record_database_versions.sh"
MANIFEST = ROOT / "conf/databases.tsv"

sys.path.insert(0, str(ROOT / "scripts"))
import finalize_provenance as fp  # noqa: E402

# Every database Ludmil named, plus the host-depletion indexes.
EXPECTED_IN_MANIFEST = ("kraken_db", "genomad_db", "gtdbtk_db", "checkm2_db",
                        "hg38_db", "t2t_phix_db")

MANIFEST_COLUMNS = 12


def read_manifest(path=MANIFEST):
    rows = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        rows[fields[0]] = fields
    return rows


def run_helper(out, *pairs, manifest=None):
    command = ["bash", str(HELPER), str(out)]
    if manifest:
        command += ["--manifest", str(manifest)]
    result = subprocess.run(command + list(pairs), capture_output=True, text=True)
    record = {}
    if Path(out).exists():
        with Path(out).open(newline="") as handle:
            record = {row["name"]: row for row in csv.DictReader(handle, delimiter="\t")}
    return result, record


class TheFrozenManifest(unittest.TestCase):
    def test_it_exists_and_covers_every_database_the_pipeline_reads(self):
        self.assertTrue(MANIFEST.is_file(), f"{MANIFEST} is missing")
        rows = read_manifest()
        for name in EXPECTED_IN_MANIFEST:
            self.assertIn(name, rows, f"conf/databases.tsv does not record {name}")

    def test_every_row_is_complete(self):
        for name, fields in read_manifest().items():
            self.assertEqual(len(fields), MANIFEST_COLUMNS,
                             f"{name}: {len(fields)} columns, expected {MANIFEST_COLUMNS}")
            release, size, count, method, fingerprint = fields[2], fields[4], fields[5], fields[6], fields[7]
            self.assertTrue(release.strip(), f"{name}: no release recorded")
            self.assertTrue(size.isdigit() and int(size) > 0, f"{name}: no size recorded")
            self.assertTrue(count.isdigit() and int(count) > 0, f"{name}: no file count recorded")
            self.assertIn(method, ("sha256:listing", "sha256:content"), f"{name}: {method}")
            self.assertRegex(fingerprint, r"^[0-9a-f]{64}$", f"{name}: fingerprint is not a sha256")

    def test_the_key_file_backstop_is_a_real_sha256_where_one_was_taken(self):
        # The listing fingerprint cannot see an in-place edit that preserves file sizes.
        # The named key file's content sha256 is what covers that, so where a key file is
        # named its checksum has to be a real one -- and for the databases Ludmil listed
        # by name, one must be named.
        rows = read_manifest()
        for name in ("kraken_db", "genomad_db", "gtdbtk_db", "checkm2_db"):
            key_file, checksum = rows[name][8], rows[name][9]
            self.assertNotEqual(key_file, "-", f"{name}: no key file named")
            self.assertRegex(checksum, r"^[0-9a-f]{64}$", f"{name}: key file checksum")

    def test_the_parameters_the_report_fingerprints_are_the_ones_the_manifest_freezes(self):
        rows = read_manifest()
        for name in rows:
            self.assertIn(name, fp.DATABASE_PARAMS,
                          f"conf/databases.tsv freezes {name}, which RUN_REPORT never looks at")


class TheRunTimeRecord(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / "database_versions.tsv"
        # A stand-in for a release directory, with the marker geNomad ships.
        self.db = self.root / "genomad_db_v1.9"
        (self.db / "sub").mkdir(parents=True)
        (self.db / "version.txt").write_text("1.9\n")
        (self.db / "sub" / "markers.bin").write_bytes(b"x" * 100)

    def frozen_manifest(self, fingerprint, name="genomad_db"):
        path = self.root / "databases.tsv"
        path.write_text("# name\tdatabase\trelease\tpath\tsize\tcount\tmethod\tfingerprint\n"
                        + "\t".join([name, "geNomad", "1.9", str(self.db),
                                     "100", "2", "sha256:listing", fingerprint]) + "\n")
        return path

    def test_it_reads_the_release_the_database_declares_not_the_directory_name(self):
        _result, record = run_helper(self.out, f"genomad_db={self.db}")
        row = record["genomad_db"]
        self.assertEqual(row["version"], "1.9")
        self.assertEqual(row["status"], "ok")
        self.assertEqual(row["kind"], "directory")
        self.assertEqual(row["file_count"], "2")
        self.assertEqual(row["fingerprint_method"], "sha256:listing")
        self.assertRegex(row["fingerprint"], r"^[0-9a-f]{64}$")

    def test_a_database_matching_the_frozen_manifest_is_reported_frozen(self):
        _result, record = run_helper(self.out, f"genomad_db={self.db}")
        manifest = self.frozen_manifest(record["genomad_db"]["fingerprint"])
        _result, record = run_helper(self.out, f"genomad_db={self.db}", manifest=manifest)
        self.assertEqual(record["genomad_db"]["known"], "frozen")

    def test_a_rebuilt_database_is_reported_as_changed(self):
        # The whole point: the path is identical, the contents are not.
        _result, record = run_helper(self.out, f"genomad_db={self.db}")
        manifest = self.frozen_manifest(record["genomad_db"]["fingerprint"])
        before = record["genomad_db"]["fingerprint"]
        (self.db / "sub" / "markers.bin").write_bytes(b"y" * 250)   # same name, new size
        _result, record = run_helper(self.out, f"genomad_db={self.db}", manifest=manifest)
        self.assertNotEqual(record["genomad_db"]["fingerprint"], before)
        self.assertTrue(record["genomad_db"]["known"].startswith("CHANGED"),
                        record["genomad_db"]["known"])

    def test_an_added_file_changes_the_fingerprint(self):
        _result, record = run_helper(self.out, f"genomad_db={self.db}")
        before = record["genomad_db"]["fingerprint"]
        (self.db / "extra.bin").write_bytes(b"z")
        _result, record = run_helper(self.out, f"genomad_db={self.db}")
        self.assertNotEqual(record["genomad_db"]["fingerprint"], before)

    def test_mtime_alone_does_not_change_the_fingerprint(self):
        # An rsync or a restore from backup rewrites mtimes; that is the same database.
        _result, record = run_helper(self.out, f"genomad_db={self.db}")
        before = record["genomad_db"]["fingerprint"]
        import os
        os.utime(self.db / "sub" / "markers.bin", (0, 0))
        _result, record = run_helper(self.out, f"genomad_db={self.db}")
        self.assertEqual(record["genomad_db"]["fingerprint"], before)

    def test_a_database_the_manifest_has_never_seen(self):
        manifest = self.frozen_manifest("0" * 64, name="something_else")
        _result, record = run_helper(self.out, f"genomad_db={self.db}", manifest=manifest)
        self.assertEqual(record["genomad_db"]["known"], "not_in_manifest")

    def test_a_small_file_gets_a_real_content_checksum(self):
        small = self.root / "index.dmnd"
        small.write_bytes(b"abc")
        _result, record = run_helper(self.out, f"checkm2_db={small}")
        self.assertEqual(record["checkm2_db"]["fingerprint_method"], "sha256:content")
        # sha256("abc")
        self.assertEqual(record["checkm2_db"]["fingerprint"],
                         "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")

    def test_an_unset_parameter_is_skipped_not_recorded_as_missing(self):
        # The caller passes every database unconditionally; most runs set few of them.
        result, record = run_helper(self.out, "kraken_db=", "gtdbtk_db=null",
                                    f"genomad_db={self.db}")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(set(record), {"genomad_db"})

    def test_a_path_that_is_not_there_is_recorded_not_fatal(self):
        result, record = run_helper(self.out, f"gtdbtk_db={self.root}/nope")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record["gtdbtk_db"]["status"], "missing")

    def test_it_never_fails_the_run(self):
        result, _record = run_helper(self.out, f"gtdbtk_db={self.root}/nope",
                                     "kraken_db=/root/forbidden")
        self.assertEqual(result.returncode, 0, result.stderr)


class TheRunReport(unittest.TestCase):
    def rows(self, **overrides):
        row = {"name": "gtdbtk_db", "path": "/dbs/gtdbtk_r232", "kind": "directory",
               "status": "ok", "version": "r232", "size_bytes": "100482230343",
               "file_count": "264", "fingerprint_method": "sha256:listing",
               "fingerprint": "e" * 64, "known": "frozen"}
        row.update(overrides)
        return {"databases": [row]}

    def test_it_prints_the_release_and_the_path(self):
        text = "\n".join(fp.render_database_versions(self.rows()))
        self.assertIn("DATABASE VERSIONS (as read)", text)
        self.assertIn("r232", text)
        self.assertIn("/dbs/gtdbtk_r232", text)
        self.assertIn("100.5 GB", text)
        self.assertIn("frozen", text)

    def test_a_changed_database_is_called_out_in_words(self):
        text = "\n".join(fp.render_database_versions(self.rows(known="CHANGED(frozen abc...)")))
        self.assertIn("DATABASE CHANGED", text)
        self.assertIn("were NOT produced against the frozen database build", text)

    def test_an_unreadable_database_is_called_out(self):
        text = "\n".join(fp.render_database_versions(self.rows(status="unreadable")))
        self.assertIn("not readable at finalisation", text)

    def test_a_run_with_no_database_says_so_rather_than_printing_nothing(self):
        text = "\n".join(fp.render_database_versions({}))
        self.assertIn("not pointed at any reference database", text)

    def test_a_failure_to_record_is_reported_not_swallowed(self):
        text = "\n".join(fp.render_database_versions({"error": "boom"}))
        self.assertIn("boom", text)


class TheCollector(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.db = self.root / "genomad_db_v1.9"
        self.db.mkdir()
        (self.db / "version.txt").write_text("1.9\n")

    def test_it_writes_the_record_next_to_the_other_run_artefacts(self):
        summary = fp.collect_database_versions({"genomad_db": str(self.db)},
                                               ROOT, self.root / "runs")
        self.assertEqual([row["name"] for row in summary["databases"]], ["genomad_db"])
        self.assertTrue((self.root / "runs" / "database_versions.tsv").is_file())
        # It is compared against the repo's own frozen manifest -- and this stand-in
        # directory is not the frozen geNomad 1.9, so the comparison has to say so.
        self.assertEqual(summary["manifest"], str(MANIFEST))
        self.assertTrue(summary["databases"][0]["known"].startswith("CHANGED"),
                        summary["databases"][0]["known"])

    def test_a_run_with_no_database_parameters_collects_nothing(self):
        self.assertEqual(fp.collect_database_versions({"sample_type": "wgs"},
                                                      ROOT, self.root / "runs"), {})

    def test_the_summary_survives_a_round_trip_through_the_json_record(self):
        summary = fp.collect_database_versions({"genomad_db": str(self.db)},
                                               ROOT, self.root / "runs")
        self.assertEqual(json.loads(json.dumps(summary)), summary)


if __name__ == "__main__":
    unittest.main()
