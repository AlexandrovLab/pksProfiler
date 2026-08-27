#!/usr/bin/env python3

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
CLB_GENES = tuple(f"clb{letter}" for letter in "ABCDEFGHIJKLMNOPQRS")


def load_script(name):
    path = REPO / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReferenceIntegrityTests(unittest.TestCase):
    def test_gff_contains_exactly_clb_a_through_clb_s(self):
        annotation = REPO / "ref" / "annotations" / "IHE3034.clbA-clbS.gff"
        genes = []

        with annotation.open() as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "gene":
                    continue

                attributes = dict(
                    field.split("=", 1)
                    for field in fields[8].split(";")
                    if "=" in field
                )
                name = attributes.get("Name")
                if name in CLB_GENES:
                    genes.append(name)

        self.assertEqual(sorted(genes), sorted(CLB_GENES))
        self.assertEqual(len(genes), 19)

    def test_hmm_database_contains_exactly_19_clb_models(self):
        hmm = REPO / "ref" / "hmm" / "clb_all_dna.hmm"
        models = []

        with hmm.open() as handle:
            for line in handle:
                if line.startswith("NAME"):
                    model = line.split()[1].removesuffix(".cds.aln")
                    models.append(model)

        self.assertEqual(sorted(models), sorted(CLB_GENES))
        self.assertEqual(len(models), 19)

    def test_complete_bowtie2_index_is_present(self):
        prefix = (
            REPO
            / "indices"
            / "GCF_000025745.1"
            / "GCF_000025745.1_ASM2574v1_genomic"
        )
        suffixes = ("1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2")

        for suffix in suffixes:
            index_file = Path(f"{prefix}.{suffix}")
            self.assertTrue(index_file.is_file(), index_file)
            self.assertGreater(index_file.stat().st_size, 0, index_file)


class TaxonomyJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.matrix_module = load_script("build_clb_species_matrix")

    def test_mate_identifiers_remain_distinct(self):
        self.assertEqual(self.matrix_module.normalize_read_id("pair/1"), "pair/1")
        self.assertEqual(self.matrix_module.normalize_read_id("pair/2"), "pair/2")

    def test_one_mates_species_is_not_assigned_to_the_other_mate(self):
        with tempfile.TemporaryDirectory() as directory:
            read_gene = Path(directory) / "read_gene.tsv"
            read_gene.write_text(
                "read_id\tGene\toverlap_bp\n"
                "pair/1\tclbA\t100\n"
                "pair/2\tclbB\t100\n"
            )

            matrix = self.matrix_module.build_matrix(
                read_gene,
                {"pair/1": {"562"}},
                {"pair/1"},
                {"562": "Escherichia coli"},
            )

        self.assertEqual(matrix[("Escherichia coli", "562")]["clbA"], 1)
        self.assertEqual(matrix[("Escherichia coli", "562")]["clbB"], 0)
        self.assertEqual(matrix[("Unclassified", "0")]["clbB"], 1)


class CombinedSpeciesSupportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.combine_module = load_script("combine_clb_species_support")

    def test_combined_output_keeps_sample_and_19_gene_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            input_path = directory / "sample1.clb_species_support.tsv"
            output_path = directory / "combined.tsv"
            gene_values = [1, *([0] * 18)]

            with input_path.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["Species", "TaxID", *CLB_GENES, "Total"])
                writer.writerow(["Escherichia coli", "562", *gene_values, 1])

            rows = self.combine_module.load_rows([input_path])
            self.combine_module.write_output(rows, output_path)

            with output_path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                combined = list(reader)

        self.assertEqual(tuple(reader.fieldnames), ("Sample", "Species", "TaxID", *CLB_GENES, "Total"))
        self.assertEqual(combined[0]["Sample"], "sample1")
        self.assertEqual(combined[0]["clbA"], "1")
        self.assertEqual(combined[0]["Total"], "1")


class QCSummaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qc_module = load_script("build_qc_summary")

    @staticmethod
    def write_fragment(path, sample, rows):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Sample", "Metric", "Value"])
            for metric, value in rows:
                writer.writerow([sample, metric, value])

    def test_bam_qc_summary_reports_all_alignment_stages(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            sample = "bam_sample"
            fragments = []
            stages = [
                [("bam_input_primary_records", 1000), ("extracted_unmapped_reads", 200)],
                [("filter_input_reads", 200), ("reads_after_fastp", 180)],
                [("reads_after_hg38", 120), ("reads_after_t2t_phix", 110)],
                [
                    ("reads_mapping_clb", 18),
                    ("clb_genes_detected", 12),
                ],
            ]

            for index, rows in enumerate(stages):
                fragment = directory / f"stage{index}.tsv"
                self.write_fragment(fragment, sample, rows)
                fragments.append(fragment)

            output = directory / "summary.tsv"
            metrics = self.qc_module.load_fragments(fragments)
            self.qc_module.write_summary(metrics, output)

            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(row["input_reads"], "1000")
        self.assertEqual(row["unmapped_reads"], "200")
        self.assertEqual(row["reads_after_fastp"], "180")
        self.assertEqual(row["reads_after_hg38"], "120")
        self.assertEqual(row["reads_after_t2t_phix"], "110")
        self.assertEqual(row["reads_mapping_clb"], "18")
        self.assertEqual(row["clb_genes_detected"], "12")
        self.assertEqual(row["reads_passing_hmm_threshold"], "NA")
        self.assertEqual(row["clb_genes_detected_hmm"], "NA")

    def test_fastq_hmm_only_summary_uses_input_count_and_na_alignment_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            sample = "fastq_sample"
            filter_fragment = directory / "filter.tsv"
            depletion_fragment = directory / "depletion.tsv"
            hmm_fragment = directory / "hmm.tsv"
            output = directory / "summary.tsv"

            self.write_fragment(
                filter_fragment,
                sample,
                [("filter_input_reads", 400), ("reads_after_fastp", 360)],
            )
            self.write_fragment(
                depletion_fragment,
                sample,
                [("reads_after_hg38", 250), ("reads_after_t2t_phix", 230)],
            )
            self.write_fragment(
                hmm_fragment,
                sample,
                [
                    ("reads_passing_hmm_threshold", 17),
                    ("clb_genes_detected_hmm", 8),
                ],
            )

            metrics = self.qc_module.load_fragments(
                [filter_fragment, depletion_fragment, hmm_fragment]
            )
            self.qc_module.write_summary(metrics, output)

            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(row["input_reads"], "400")
        self.assertEqual(row["unmapped_reads"], "400")
        self.assertEqual(row["reads_mapping_clb"], "NA")
        self.assertEqual(row["clb_genes_detected"], "NA")
        self.assertEqual(row["reads_passing_hmm_threshold"], "17")
        self.assertEqual(row["clb_genes_detected_hmm"], "8")

    def test_qc_summary_rejects_increasing_downstream_read_count(self):
        values = {
            "filter_input_reads": 100,
            "reads_after_fastp": 90,
            "reads_after_hg38": 95,
            "reads_after_t2t_phix": 80,
        }

        with self.assertRaisesRegex(ValueError, "Impossible QC counts"):
            self.qc_module.output_row("invalid_sample", values)


class ExampleOutputTests(unittest.TestCase):
    def test_all_example_count_tables_have_19_ordered_integer_rows(self):
        example_files = sorted((REPO / "examples" / "results").glob("*/*.txt"))
        self.assertEqual(len(example_files), 4)

        for path in example_files:
            with self.subTest(path=path):
                with path.open(newline="") as handle:
                    rows = list(csv.reader(handle, delimiter="\t"))

                self.assertEqual(rows[0][0], "Gene")
                self.assertEqual(tuple(row[0] for row in rows[1:]), CLB_GENES)
                self.assertEqual(len(rows), 20)

                for row in rows[1:]:
                    self.assertGreaterEqual(int(row[1]), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
