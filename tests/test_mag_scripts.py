#!/usr/bin/env python3
"""Tests for MAG-module helper scripts."""

import csv
import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = str(REPO / "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


def load_script(name):
    path = REPO / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# extract_genomic_context
# ---------------------------------------------------------------------------

GFF_CONTENT = """\
##gff-version 3
##sequence-region contig1 1 200000
contig1\tProdigal\tCDS\t1000\t2000\t.\t+\t.\tID=PROKKA_00001;gene=intA;product=integrase
contig1\tProdigal\tCDS\t3000\t4000\t.\t+\t.\tID=PROKKA_00002;gene=clbA_prot;product=colibactin protein A
contig1\tProdigal\tCDS\t100000\t101000\t.\t-\t.\tID=PROKKA_00003;gene=tnpA;product=transposase
contig2\tProdigal\tCDS\t500\t1500\t.\t+\t.\tID=PROKKA_00004;gene=gyrB;product=gyrase subunit B
"""

GFF_WITH_GENE_FEATURES = """\
##gff-version 3
##sequence-region contig1 1 200000
contig1\tProkka\tgene\t1000\t2000\t.\t+\t.\tID=gene_1
contig1\tProkka\tCDS\t1000\t2000\t.\t+\t.\tID=PROKKA_00001;gene=intA;product=integrase
contig1\tProkka\tgene\t3000\t4000\t.\t+\t.\tID=gene_2
contig1\tProkka\tCDS\t3000\t4000\t.\t+\t.\tID=PROKKA_00002;gene=clbA;product=colibactin A
"""

TBLOUT_CONTENT = """\
#                                                               --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias  exp  dom  seq  from    to  from    to  from    to  acc description of target
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ ----- ---- ---- ---- ----- ----- ----- ----- ----- ----- ---- ---------------------
PROKKA_00002         -          clbA                 -            1.2e-10  35.0   0.0   1.5e-10   34.8   0.0   1.0    1    1     1   200    10   205    8   207  0.95 colibactin A
"""


class TestExtractGenomicContext(unittest.TestCase):
    def setUp(self):
        self.mod = load_script("extract_genomic_context")

    def _write(self, tmp, name, content):
        p = os.path.join(tmp, name)
        with open(p, "w") as f:
            f.write(content)
        return p

    def test_parse_prokka_gff_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = self._write(tmp, "test.gff", GFF_CONTENT)
            genes = self.mod.parse_prokka_gff(gff)
        self.assertEqual(len(genes), 4)

    def test_parse_prokka_gff_excludes_gene_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = self._write(tmp, "test.gff", GFF_WITH_GENE_FEATURES)
            genes = self.mod.parse_prokka_gff(gff)
        self.assertEqual(len(genes), 2)

    def test_parse_prokka_gff_gene_names_from_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = self._write(tmp, "test.gff", GFF_WITH_GENE_FEATURES)
            genes = self.mod.parse_prokka_gff(gff)
        self.assertEqual(genes[0]["locus_tag"], "PROKKA_00001")
        self.assertEqual(genes[0]["gene"], "intA")

    def test_parse_hmmsearch_tblout(self):
        with tempfile.TemporaryDirectory() as tmp:
            tbl = self._write(tmp, "test.tblout", TBLOUT_CONTENT)
            hits = self.mod.parse_hmmsearch_tblout(tbl, evalue_threshold=1e-5)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["clb_gene"], "clbA")

    def test_is_integrase(self):
        self.assertTrue(self.mod.is_integrase({"gene": "intA", "product": "integrase"}))
        self.assertFalse(self.mod.is_integrase({"gene": "gyrB", "product": "gyrase"}))

    def test_is_transposase(self):
        self.assertTrue(self.mod.is_transposase({"gene": "tnpA", "product": "transposase"}))
        self.assertFalse(self.mod.is_transposase({"gene": "gyrB", "product": "gyrase"}))

    def test_get_flanking_within_window(self):
        genes = [
            {"contig": "c1", "start": 1000, "end": 2000, "locus_tag": "L1", "gene": "g1", "product": ""},
            {"contig": "c1", "start": 80000, "end": 81000, "locus_tag": "L2", "gene": "g2", "product": ""},
            {"contig": "c2", "start": 1000, "end": 2000, "locus_tag": "L3", "gene": "g3", "product": ""},
        ]
        result = self.mod.get_flanking(genes, "c1", 5000, window=50000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["locus_tag"], "L1")

    def test_end_to_end_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = self._write(tmp, "bin.gff", GFF_CONTENT)
            tbl = self._write(tmp, "bin.tblout", TBLOUT_CONTENT)
            out = os.path.join(tmp, "out.tsv")
            import sys
            sys.argv = ["extract_genomic_context.py",
                        "--gff", gff, "--tblout", tbl,
                        "--evalue", "1e-5", "--out", out]
            self.mod.main()
            with open(out) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["clb_gene"], "clbA")
        # integrase (PROKKA_00001) is within 50kb of PROKKA_00002 on contig1
        self.assertEqual(rows[0]["has_integrase"], "True")
        self.assertEqual(rows[0]["has_transposase"], "False")


# ---------------------------------------------------------------------------
# build_mag_summary
# ---------------------------------------------------------------------------

# Column names exactly as CheckM2 writes them (checkm2/predictQuality.py). The
# previous fixture said "Genome_size", matching the bug in the reader rather than the
# tool -- which is how M4 stayed hidden while this test passed.
CHECKM2_CONTENT = """\
Name\tCompleteness\tContamination\tGenome_Size\tContig_N50
bin_001\t92.5\t1.2\t4500000\t125000
bin_002\t45.0\t8.0\t2100000\t18000
"""

GTDBTK_CONTENT = """\
user_genome\tclassification
bin_001\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
bin_002\td__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus acidophilus
"""

TBLOUT_PKS_CONTENT = """\
#
bin_001_1\t-\tclbA\t-\t1e-10\t35.0\t0.0\t1.5e-10\t34.8\t0.0\t1.0\t1\t1\t1\t200\t10\t205\t8\t207\t0.95\thit
bin_001_2\t-\tclbB\t-\t2e-08\t30.0\t0.0\t2.5e-08\t29.8\t0.0\t1.0\t1\t1\t1\t200\t10\t205\t8\t207\t0.95\thit
"""

CONTEXT_CONTENT = """\
locus_tag\tclb_gene\tevalue\tcontig\thas_integrase\thas_transposase\tflanking_genes
bin_001_1\tclbA\t1e-10\tbin_001\tTrue\tFalse\tgyrB;fliA
"""

# magBinLocusEvidence's output for a bin whose own assembly aligns across the
# canonical locus -- the alignment-confirmed call build_mag_summary.py and
# build_community_prophage.py now join in, instead of trusting HMM gene count alone.
LOCUS_EVIDENCE_POSITIVE = """\
sample\tbin_id\tlocus_tier\tlocus_genes_detected\tlocus_genes\tlocus_breadth
SAMPLE1\tbin_001\tmulti_gene\t5\tclbA,clbB,clbD,clbP,clbQ\t0.020000
"""


class TestBuildMagSummary(unittest.TestCase):
    def setUp(self):
        self.mod = load_script("build_mag_summary")

    def _write(self, tmp, name, content):
        p = os.path.join(tmp, name)
        with open(p, "w") as f:
            f.write(content)
        return p

    def test_parse_checkm2(self):
        with tempfile.TemporaryDirectory() as tmp:
            tsv = self._write(tmp, "quality.tsv", CHECKM2_CONTENT)
            result = self.mod.parse_checkm2(tsv)
        self.assertAlmostEqual(result["bin_001"]["completeness"], 92.5)
        self.assertAlmostEqual(result["bin_001"]["contamination"], 1.2)

    def test_parse_gtdbtk(self):
        with tempfile.TemporaryDirectory() as tmp:
            tsv = self._write(tmp, "gtdbtk.tsv", GTDBTK_CONTENT)
            result = self.mod.parse_gtdbtk(tsv)
        self.assertIn("Enterobacterales", result["bin_001"])
        self.assertNotIn("Enterobacterales", result["bin_002"])

    def test_unexpected_taxon_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            gtdbtk = self._write(tmp, "gtdbtk.tsv", GTDBTK_CONTENT)
            result = self.mod.parse_gtdbtk(gtdbtk)
        self.assertFalse(self.mod.is_unexpected_taxon(result["bin_001"]))
        self.assertTrue(self.mod.is_unexpected_taxon(result["bin_002"]))

    def test_is_unexpected_taxon_unclassified(self):
        self.assertFalse(self.mod.is_unexpected_taxon("unclassified"))

    def test_is_unexpected_taxon_classified_other_order(self):
        self.assertTrue(self.mod.is_unexpected_taxon(
            "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;"
            "f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus acidophilus"
        ))

    def test_is_unexpected_taxon_enterobacterales(self):
        self.assertFalse(self.mod.is_unexpected_taxon(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"
        ))

    def test_is_unexpected_taxon_empty_string(self):
        self.assertFalse(self.mod.is_unexpected_taxon(""))

    def test_end_to_end_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            checkm2 = self._write(tmp, "quality.tsv", CHECKM2_CONTENT)
            gtdbtk = self._write(tmp, "gtdbtk.tsv", GTDBTK_CONTENT)
            tblout_dir = os.path.join(tmp, "tblout")
            os.makedirs(tblout_dir)
            self._write(tblout_dir, "bin_001.tblout", TBLOUT_PKS_CONTENT)
            context_dir = os.path.join(tmp, "context")
            os.makedirs(context_dir)
            self._write(context_dir, "bin_001.context.tsv", CONTEXT_CONTENT)
            locus_dir = os.path.join(tmp, "locus")
            os.makedirs(locus_dir)
            self._write(locus_dir, "bin_001.locus_evidence.tsv", LOCUS_EVIDENCE_POSITIVE)
            out = os.path.join(tmp, "summary.tsv")
            import sys
            sys.argv = ["build_mag_summary.py",
                        "--checkm2", checkm2, "--gtdbtk", gtdbtk,
                        "--tblout_dir", tblout_dir, "--context_dir", context_dir,
                        "--locus_dir", locus_dir,
                        "--sample", "SAMPLE1", "--evalue", "1e-5", "--out", out]
            self.mod.main()
            with open(out) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sample"], "SAMPLE1")
        self.assertEqual(rows[0]["bin_id"], "bin_001")
        self.assertEqual(rows[0]["clb_genes_detected"], "2")
        self.assertIn("clbA", rows[0]["clb_genes"])
        self.assertIn("clbB", rows[0]["clb_genes"])
        self.assertEqual(rows[0]["has_integrase"], "True")
        self.assertEqual(rows[0]["unexpected_taxon_flag"], "False")
        self.assertEqual(rows[0]["locus_tier"], "multi_gene")
        self.assertEqual(rows[0]["locus_genes_detected"], "5")

    def test_missing_locus_evidence_reads_as_unmeasured_not_negative(self):
        """A bin absent from locus_dir was never aligned -- NA, not a failed tier."""
        with tempfile.TemporaryDirectory() as tmp:
            checkm2 = self._write(tmp, "quality.tsv", CHECKM2_CONTENT)
            gtdbtk = self._write(tmp, "gtdbtk.tsv", GTDBTK_CONTENT)
            tblout_dir = os.path.join(tmp, "tblout")
            os.makedirs(tblout_dir)
            self._write(tblout_dir, "bin_001.tblout", TBLOUT_PKS_CONTENT)
            context_dir = os.path.join(tmp, "context")
            os.makedirs(context_dir)
            locus_dir = os.path.join(tmp, "locus")
            os.makedirs(locus_dir)
            out = os.path.join(tmp, "summary.tsv")
            import sys
            sys.argv = ["build_mag_summary.py",
                        "--checkm2", checkm2, "--gtdbtk", gtdbtk,
                        "--tblout_dir", tblout_dir, "--context_dir", context_dir,
                        "--locus_dir", locus_dir,
                        "--sample", "SAMPLE1", "--evalue", "1e-5", "--out", out]
            self.mod.main()
            with open(out) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(rows[0]["locus_tier"], "NA")

    # -----------------------------------------------------------------------
    # U1: unbinned-contig pks positivity (real bin binning left something out)
    # -----------------------------------------------------------------------

    # magBinLocusEvidence's output for the per-sample pool of contigs MetaBAT2 never
    # placed in any bin -- "unbinned" is a stand-in bin_id, exactly like a real one,
    # produced by the identical summarize_mag_bin_locus_evidence.py --bin-id call.
    LOCUS_EVIDENCE_UNBINNED_POSITIVE = (
        "sample\tbin_id\tlocus_tier\tlocus_genes_detected\tlocus_genes\tlocus_breadth\n"
        "SAMPLE1\tunbinned\tbroad_island\t9\tclbA,clbB,clbD,clbK,clbN,clbO,clbP,clbQ,clbS\t0.180000\n"
    )

    def _write_bin_and_unbinned_fixture(self, tmp, unbinned_content):
        checkm2 = self._write(tmp, "quality.tsv", CHECKM2_CONTENT)
        gtdbtk = self._write(tmp, "gtdbtk.tsv", GTDBTK_CONTENT)
        tblout_dir = os.path.join(tmp, "tblout")
        os.makedirs(tblout_dir)
        self._write(tblout_dir, "bin_001.tblout", TBLOUT_PKS_CONTENT)
        context_dir = os.path.join(tmp, "context")
        os.makedirs(context_dir)
        locus_dir = os.path.join(tmp, "locus")
        os.makedirs(locus_dir)
        self._write(locus_dir, "bin_001.locus_evidence.tsv", LOCUS_EVIDENCE_POSITIVE)
        if unbinned_content is not None:
            self._write(locus_dir, "unbinned.locus_evidence.tsv", unbinned_content)
        return checkm2, gtdbtk, tblout_dir, context_dir, locus_dir

    def _run_with(self, tmp, checkm2, gtdbtk, tblout_dir, context_dir, locus_dir):
        out = os.path.join(tmp, "summary.tsv")
        import sys
        sys.argv = ["build_mag_summary.py",
                    "--checkm2", checkm2, "--gtdbtk", gtdbtk,
                    "--tblout_dir", tblout_dir, "--context_dir", context_dir,
                    "--locus_dir", locus_dir,
                    "--sample", "SAMPLE1", "--evalue", "1e-5", "--out", out]
        self.mod.main()
        with open(out) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_unbinned_positive_call_is_its_own_row_not_folded_into_a_bin(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = self._write_bin_and_unbinned_fixture(tmp, self.LOCUS_EVIDENCE_UNBINNED_POSITIVE)
            rows = self._run_with(tmp, *fixture)
        self.assertEqual(len(rows), 2)
        bins = [r for r in rows if r["unit_type"] == "bin"]
        unbinned = [r for r in rows if r["unit_type"] == "unbinned"]
        self.assertEqual(len(bins), 1)
        self.assertEqual(bins[0]["bin_id"], "bin_001")
        self.assertEqual(len(unbinned), 1)
        self.assertEqual(unbinned[0]["bin_id"], "unbinned")
        self.assertEqual(unbinned[0]["locus_tier"], "broad_island")
        self.assertEqual(unbinned[0]["locus_genes_detected"], "9")
        # Not a genome: none of the bin-specific fields were measured for it.
        self.assertEqual(unbinned[0]["taxonomy"], "NA")
        self.assertEqual(unbinned[0]["completeness"], "NA")
        self.assertEqual(unbinned[0]["clb_genes_detected"], "NA")

    def test_existing_bin_rows_are_explicitly_typed_bin(self):
        """Old readers keyed only on bin_id; unit_type must not silently change
        what a plain per-bin row looks like."""
        with tempfile.TemporaryDirectory() as tmp:
            fixture = self._write_bin_and_unbinned_fixture(tmp, None)
            rows = self._run_with(tmp, *fixture)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["unit_type"], "bin")

    def test_no_unbinned_evidence_file_means_no_unbinned_row(self):
        """Nothing to report when metabat2Bin.out.unbinned never emitted for this
        sample (e.g. every contig got binned, or there were no contigs at all)."""
        with tempfile.TemporaryDirectory() as tmp:
            fixture = self._write_bin_and_unbinned_fixture(tmp, None)
            rows = self._run_with(tmp, *fixture)
        self.assertNotIn("unbinned", [r["bin_id"] for r in rows])

    def test_a_negative_unbinned_tier_is_still_reported(self):
        """Checked and found nothing is not the same as never checked -- the row
        still appears, distinguishing the two the same way magSampleStatus does."""
        negative = (
            "sample\tbin_id\tlocus_tier\tlocus_genes_detected\tlocus_genes\tlocus_breadth\n"
            "SAMPLE1\tunbinned\tnegative\t0\t\t0.000000\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            fixture = self._write_bin_and_unbinned_fixture(tmp, negative)
            rows = self._run_with(tmp, *fixture)
        unbinned = next(r for r in rows if r["unit_type"] == "unbinned")
        self.assertEqual(unbinned["locus_tier"], "negative")


# ---------------------------------------------------------------------------
# build_community_prophage
# ---------------------------------------------------------------------------

# A realistic provirus: the 401 bp interval this fixture used to carry is below the
# T3 noise floor, and the length was incidental to what these tests check (that only
# `provirus` topology is an inventory candidate). The floor itself is covered in
# tests/test_provirus_noise_floor.py.
GENOMAD_CONTENT = """\
seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy
contig1|provirus_100_42100\t42001\tProvirus\t100-42100\t48\t11\t0.97\tNA\t3\t10\tViruses;Caudoviricetes
free_virus\t35000\tNo terminal repeats\tNA\t40\t11\t0.99\tNA\t4\t12\tViruses
"""


class TestCommunityProphage(unittest.TestCase):
    def setUp(self):
        self.mod = load_script("build_community_prophage")

    def _write(self, directory, name, content):
        path = os.path.join(directory, name)
        with open(path, "w") as fh:
            fh.write(content)
        return path

    def test_only_integrated_proviruses_are_inventory_candidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = self._write(tmp, "bin_001.virus_summary.tsv", GENOMAD_CONTENT)
            rows = self.mod.read_prophages(summary)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["host_contig"], "contig1")
        self.assertEqual((rows[0]["start"], rows[0]["end"]), (100, 42100))

    def test_sos_markers_come_from_explicit_gene_names(self):
        gff = """##gff-version 3
c1\tProkka\tCDS\t1\t100\t.\t+\t.\tID=x1;gene=recA;product=recombinase A
c1\tProkka\tCDS\t200\t300\t.\t+\t.\tID=x2;gene=lexA;product=LexA repressor
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(tmp, "bin.gff", gff)
            self.assertEqual(self.mod.annotation_flags(path), (True, True))

    def test_interaction_table_is_hypothesis_not_induction_call(self):
        gff = """##gff-version 3
c1\tProkka\tCDS\t1\t100\t.\t+\t.\tID=x1;gene=recA;product=recombinase A
c1\tProkka\tCDS\t200\t300\t.\t+\t.\tID=x2;gene=lexA;product=LexA repressor
"""
        three_hits = """#
x1\t-\tclbB\t-\t1e-20\t50\t0\t1e-20\t50\t0\t1\t1\t1\t1\t10\t1\t10\t1\t10\t1\thit
x2\t-\tclbC\t-\t1e-20\t50\t0\t1e-20\t50\t0\t1\t1\t1\t1\t10\t1\t10\t1\t10\t1\thit
x3\t-\tclbS\t-\t1e-20\t50\t0\t1e-20\t50\t0\t1\t1\t1\t1\t10\t1\t10\t1\t10\t1\thit
"""
        with tempfile.TemporaryDirectory() as tmp:
            gff_dir = os.path.join(tmp, "gff")
            tbl_dir = os.path.join(tmp, "tbl")
            virus_dir = os.path.join(tmp, "virus")
            locus_dir = os.path.join(tmp, "locus")
            os.makedirs(gff_dir); os.makedirs(tbl_dir); os.makedirs(virus_dir); os.makedirs(locus_dir)
            self._write(gff_dir, "bin_001.gff", gff)
            self._write(tbl_dir, "bin_001.tblout", three_hits)
            self._write(virus_dir, "bin_001.virus_summary.tsv", GENOMAD_CONTENT)
            # M1: a producer now also needs its own assembly to align across the
            # canonical locus at a positive tier -- clbB/clbC alone would not qualify.
            self._write(locus_dir, "bin_001.locus_evidence.tsv", LOCUS_EVIDENCE_POSITIVE)
            taxonomy = self._write(
                tmp, "taxonomy.tsv",
                "user_genome\tclassification\nbin_001\td__Bacteria;o__Enterobacterales\n",
            )
            inventory = os.path.join(tmp, "inventory.tsv")
            interactions = os.path.join(tmp, "interactions.tsv")
            mobility = os.path.join(tmp, "mobility.tsv")
            old_argv = sys.argv
            try:
                sys.argv = [
                    "build_community_prophage.py", "--sample", "S1",
                    "--gtdbtk", taxonomy, "--gff-dir", gff_dir,
                    "--tblout-dir", tbl_dir, "--genomad-dir", virus_dir,
                    "--locus-dir", locus_dir,
                    "--prophage-out", inventory, "--interaction-out", interactions,
                    "--mobility-out", mobility,
                ]
                self.mod.main()
            finally:
                sys.argv = old_argv
            with open(interactions) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            with open(mobility) as fh:
                mobility_rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["recipient_clbS_like"], "True")
        self.assertEqual(rows[0]["susceptibility_hypothesis"], "lower")
        self.assertEqual(rows[0]["interpretation"], "co-occurrence_only_not_evidence_of_induction")
        self.assertEqual(mobility_rows[0]["provisional_pks_class"], "partial_pks_candidate")
        self.assertEqual(mobility_rows[0]["mobility_interpretation"], "no_local_HGT_marker_detected")
        self.assertEqual(mobility_rows[0]["locus_tier"], "multi_gene")

    def test_megasynthase_only_bin_is_not_a_producer(self):
        """M1: clbB/clbK/clbH clear BIOSYNTHETIC_CLB and the gene-count floor, but the
        bin's own assembly does not align across the canonical locus (the ERR525841
        Bifidobacterium pattern) -- must not produce an interaction row, and must read
        as HMM-only, not a confirmed pks candidate.
        """
        gff = """##gff-version 3
c1\tProkka\tCDS\t1\t100\t.\t+\t.\tID=x1;gene=hypothetical;product=hypothetical protein
"""
        megasynthase_only = """#
x1\t-\tclbB\t-\t1e-20\t50\t0\t1e-20\t50\t0\t1\t1\t1\t1\t10\t1\t10\t1\t10\t1\thit
x2\t-\tclbK\t-\t1e-20\t50\t0\t1e-20\t50\t0\t1\t1\t1\t1\t10\t1\t10\t1\t10\t1\thit
x3\t-\tclbH\t-\t1e-20\t50\t0\t1e-20\t50\t0\t1\t1\t1\t1\t10\t1\t10\t1\t10\t1\thit
"""
        locus_negative = (
            "sample\tbin_id\tlocus_tier\tlocus_genes_detected\tlocus_genes\tlocus_breadth\n"
            "S1\tbin_002\tnegative\t0\t\t0.000000\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            gff_dir = os.path.join(tmp, "gff")
            tbl_dir = os.path.join(tmp, "tbl")
            virus_dir = os.path.join(tmp, "virus")
            locus_dir = os.path.join(tmp, "locus")
            os.makedirs(gff_dir); os.makedirs(tbl_dir); os.makedirs(virus_dir); os.makedirs(locus_dir)
            self._write(gff_dir, "bin_002.gff", gff)
            self._write(tbl_dir, "bin_002.tblout", megasynthase_only)
            self._write(virus_dir, "bin_002.virus_summary.tsv", GENOMAD_CONTENT)
            self._write(locus_dir, "bin_002.locus_evidence.tsv", locus_negative)
            taxonomy = self._write(
                tmp, "taxonomy.tsv",
                "user_genome\tclassification\n"
                "bin_002\td__Bacteria;p__Actinobacteria;o__Bifidobacteriales;"
                "g__Bifidobacterium;s__Bifidobacterium longum\n",
            )
            inventory = os.path.join(tmp, "inventory.tsv")
            interactions = os.path.join(tmp, "interactions.tsv")
            mobility = os.path.join(tmp, "mobility.tsv")
            old_argv = sys.argv
            try:
                sys.argv = [
                    "build_community_prophage.py", "--sample", "S1",
                    "--gtdbtk", taxonomy, "--gff-dir", gff_dir,
                    "--tblout-dir", tbl_dir, "--genomad-dir", virus_dir,
                    "--locus-dir", locus_dir,
                    "--prophage-out", inventory, "--interaction-out", interactions,
                    "--mobility-out", mobility,
                ]
                self.mod.main()
            finally:
                sys.argv = old_argv
            with open(interactions) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            with open(mobility) as fh:
                mobility_rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(rows, [])  # bin_002 hits a recipient (itself) but is no producer
        self.assertEqual(mobility_rows[0]["provisional_pks_class"], "hmm_candidate_locus_not_confirmed")
        self.assertEqual(mobility_rows[0]["biosynthetic_clb_detected"], "True")
        self.assertEqual(mobility_rows[0]["locus_tier"], "negative")


# ---------------------------------------------------------------------------
# summarize_mag_bin_locus_evidence / mag_utils.read_locus_evidence
# ---------------------------------------------------------------------------

# island-start=0, island-end=1000 in these tests (real coordinates aren't the point).
LOCUS_GFF = """\
##gff-version 3
c1\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=g1;Name=clbA
c1\tRefSeq\tgene\t200\t400\t.\t+\t.\tID=g2;Name=clbB
c1\tRefSeq\tgene\t600\t900\t.\t+\t.\tID=g3;Name=clbC
"""


def paf_line(contig="c1", qlen=1000, ts=0, te=1000, aligned=1000, identity=1.0, mapq=60):
    nm = int(round(aligned * (1 - identity)))
    return "\t".join(str(x) for x in (
        contig, qlen, 0, aligned, "+", "ref", 1000, ts, te, aligned - nm, aligned, mapq,
    ))


class TestMagBinLocusEvidence(unittest.TestCase):
    def setUp(self):
        self.mod = load_script("summarize_mag_bin_locus_evidence")

    def _write(self, tmp, name, content):
        p = os.path.join(tmp, name)
        with open(p, "w") as f:
            f.write(content)
        return p

    def test_no_alignment_is_negative(self):
        with tempfile.TemporaryDirectory() as tmp:
            paf = self._write(tmp, "empty.paf", "")
            gff = self._write(tmp, "ref.gff", LOCUS_GFF)
            hits = self.mod.paf_hits(paf, 0, 1000, 0.90, 20)
        self.assertEqual(hits, [])
        self.assertEqual(self.mod.classify(0, 0.0, {
            "multi_gene": (3, .01), "broad_island": (8, .075), "extensive_island": (10, .15),
        }), "negative")

    def test_full_length_alignment_covers_all_genes_and_full_breadth(self):
        with tempfile.TemporaryDirectory() as tmp:
            paf = self._write(tmp, "full.paf", paf_line() + "\n")
            gff = self._write(tmp, "ref.gff", LOCUS_GFF)
            hits = self.mod.paf_hits(paf, 0, 1000, 0.90, 20)
            genes = self.mod.genes_covered(gff, "c1", 0, 1000, hits)
        self.assertEqual(self.mod.union_length(hits), 1000)
        self.assertEqual(genes, {"clbA", "clbB", "clbC"})

    def test_low_identity_alignment_is_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            paf = self._write(tmp, "low_identity.paf", paf_line(identity=0.5) + "\n")
            hits = self.mod.paf_hits(paf, 0, 1000, 0.90, 20)
        self.assertEqual(hits, [])

    def test_partial_alignment_covers_only_overlapping_genes(self):
        # Covers clbA (1-100) and part of clbB (200-400) but not clbC (600-900).
        with tempfile.TemporaryDirectory() as tmp:
            paf = self._write(tmp, "partial.paf", paf_line(te=300, aligned=300) + "\n")
            gff = self._write(tmp, "ref.gff", LOCUS_GFF)
            hits = self.mod.paf_hits(paf, 0, 1000, 0.90, 20)
            genes = self.mod.genes_covered(gff, "c1", 0, 1000, hits)
        self.assertEqual(genes, {"clbA", "clbB"})

    def test_classify_picks_the_highest_tier_that_holds(self):
        thresholds = {"multi_gene": (3, .01), "broad_island": (8, .075), "extensive_island": (10, .15)}
        self.assertEqual(self.mod.classify(10, 0.20, thresholds), "extensive_island")
        self.assertEqual(self.mod.classify(8, 0.08, thresholds), "broad_island")
        self.assertEqual(self.mod.classify(3, 0.02, thresholds), "multi_gene")
        self.assertEqual(self.mod.classify(1, 0.005, thresholds), "localized_indeterminate")
        self.assertEqual(self.mod.classify(0, 0.0, thresholds), "negative")

    def test_end_to_end_bifidobacterium_pattern_is_negative_not_positive(self):
        """The ERR525841 false positive: HMM hits to megasynthase domains with no
        actual alignment to the reference locus (a different genus, no shared synteny).
        """
        with tempfile.TemporaryDirectory() as tmp:
            paf = self._write(tmp, "bin.paf", "")  # no alignment survives at all
            gff = self._write(tmp, "ref.gff", LOCUS_GFF)
            out = os.path.join(tmp, "evidence.tsv")
            import sys
            sys.argv = ["summarize_mag_bin_locus_evidence.py",
                        "--sample", "ERR525841", "--bin-id", "bin.1", "--paf", paf,
                        "--gff", gff, "--contig", "c1",
                        "--island-start", "0", "--island-end", "1000", "--output", out]
            self.mod.main()
            with open(out) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(rows[0]["locus_tier"], "negative")
        self.assertEqual(rows[0]["locus_genes_detected"], "0")

    def test_read_locus_evidence_keys_by_bin_id(self):
        mag_utils = load_script("mag_utils")
        with tempfile.TemporaryDirectory() as tmp:
            self._write(tmp, "bin.1.locus_evidence.tsv",
                        "sample\tbin_id\tlocus_tier\tlocus_genes_detected\tlocus_genes\tlocus_breadth\n"
                        "S1\tbin.1\textensive_island\t12\tclbA,clbB\t0.900000\n")
            result = mag_utils.read_locus_evidence(tmp)
        self.assertIn("bin.1", result)
        self.assertEqual(result["bin.1"]["locus_tier"], "extensive_island")

    def test_bin_id_unbinned_works_unmodified_for_the_pooled_contig_pseudo_unit(self):
        """U1: Modules/pks_mag.nf reuses this exact script for the per-sample pool of
        contigs MetaBAT2 never placed in any bin, passing --bin-id unbinned instead of
        a real bin's ID. The script has no bin-specific logic -- --bin-id is an opaque
        label threaded straight through to the output's bin_id column and the output
        filename -- so this needs no code change here, only this regression test
        confirming that stays true.
        """
        with tempfile.TemporaryDirectory() as tmp:
            # A full-length alignment, same as test_full_length_alignment_covers_all_genes_
            # and_full_breadth above (100% breadth, all 3 LOCUS_GFF genes covered), but
            # scored as the unbinned pool rather than a bin. Caps at multi_gene, not a
            # higher tier, because LOCUS_GFF only has 3 genes to find -- the default
            # thresholds need >=8/>=10 genes for broad_island/extensive_island.
            paf = self._write(tmp, "unbinned.paf", paf_line() + "\n")
            gff = self._write(tmp, "ref.gff", LOCUS_GFF)
            out = os.path.join(tmp, "unbinned.locus_evidence.tsv")
            import sys
            sys.argv = ["summarize_mag_bin_locus_evidence.py",
                        "--sample", "SAMPLE1", "--bin-id", "unbinned", "--paf", paf,
                        "--gff", gff, "--contig", "c1",
                        "--island-start", "0", "--island-end", "1000", "--output", out]
            self.mod.main()
            with open(out) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(rows[0]["bin_id"], "unbinned")
            self.assertEqual(rows[0]["locus_tier"], "multi_gene")

            mag_utils = load_script("mag_utils")
            result = mag_utils.read_locus_evidence(tmp)
        self.assertIn("unbinned", result)
        self.assertEqual(result["unbinned"]["locus_tier"], "multi_gene")


if __name__ == "__main__":
    unittest.main()


class Checkm2SchemaTests(unittest.TestCase):
    """M4: a report we cannot read is a broken run, not a genome of size zero.

    build_mag_summary read `Genome_size`; CheckM2 writes `Genome_Size`. The lookup
    never matched, and `.get(..., 0)` turned that into a measurement -- 0 bp for all
    8 bins of the v0.0.2 test, which also removed the cheapest sanity check on a bin.
    """

    def setUp(self):
        self.mod = load_script("build_mag_summary")
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "quality_report.tsv"

    def test_the_real_checkm2_column_names_are_read(self):
        self.path.write_text(CHECKM2_CONTENT)
        parsed = self.mod.parse_checkm2(self.path)
        self.assertEqual(parsed["bin_001"]["genome_size"], 4500000)
        self.assertEqual(parsed["bin_001"]["contig_n50"], 125000)

    def test_a_missing_column_raises_instead_of_yielding_zero(self):
        self.path.write_text("Name\tCompleteness\tContamination\n bin_001\t92.5\t1.2\n")
        with self.assertRaises(ValueError) as caught:
            self.mod.parse_checkm2(self.path)
        self.assertIn("Genome_Size", str(caught.exception))
        self.assertIn("Columns present", str(caught.exception))

    def test_the_old_lowercase_spelling_is_not_silently_accepted(self):
        # Exactly the schema the buggy reader expected.
        self.path.write_text("Name\tCompleteness\tContamination\tGenome_size\n"
                             "bin_001\t92.5\t1.2\t4500000\n")
        with self.assertRaises(ValueError):
            self.mod.parse_checkm2(self.path)
