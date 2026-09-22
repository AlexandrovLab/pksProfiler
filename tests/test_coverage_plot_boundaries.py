"""The coverage plot keeps signal at the island edges.

F12, from Ludmil's 2026-09-19 report. plotPKS.R kept only bedgraph bins lying wholly
inside the island:

    coverage$start >= min(cytoband.df$start) & coverage$end <= max(cytoband.df$end)

A bin straddling either edge -- starting before the island and ending inside it, or
the reverse -- was dropped whole, so coverage at the first and last bins never reached
the plot. Those edges are where the mobility architecture sits (integrase at -280 bp,
tRNA-Asn at -1713 bp), so they are the least good place to lose signal.

The contig and the island offset were also written into the R while the pipeline held
them as parameters; they are arguments now.
"""
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLOT_R = ROOT / "scripts/plotPKS.R"
R_SOURCE = PLOT_R.read_text()
PLOTTING = (ROOT / "Modules/plotting.nf").read_text()

CYTOBAND = ("chrom\tstart\tend\tName\tgieStain\n"
            "NC_017628.1\t0\t1000\tclbA\tgneg\n"
            "NC_017628.1\t49000\t50767\tclbS\tgneg\n")


def code_lines(text):
    return "\n".join(l for l in text.splitlines() if not l.strip().startswith("#"))


class TheIntervalComesFromThePipeline(unittest.TestCase):
    def test_the_contig_is_not_written_into_the_script(self):
        # It was compared against a literal; now it arrives as an argument.
        self.assertNotIn('coverage$chr == \'NC_017628.1\'', code_lines(R_SOURCE))
        self.assertIn("contig        <- if (length(args) >= 4)", R_SOURCE)

    def test_the_island_offset_is_not_written_into_the_script(self):
        self.assertNotIn("+2193827", code_lines(R_SOURCE))
        self.assertIn("island_start", R_SOURCE)

    def test_the_module_passes_contig_and_interval(self):
        invocation = PLOTTING.split("plotPKS.R", 1)[1].split('"""', 1)[0]
        self.assertIn("params.pks_contig", invocation)
        self.assertIn("params.pks_shift", invocation)
        self.assertIn("params.pks_island_len", invocation)


@unittest.skipUnless(shutil.which("Rscript"), "Rscript not on PATH")
class BoundaryBinsSurvive(unittest.TestCase):
    """Runs the real script; the bin count it reports is the assertion."""

    def bins_reported(self, bedgraph_text):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d / "cyto.txt").write_text(CYTOBAND)
            (d / "cov.bedgraph").write_text(bedgraph_text)
            result = subprocess.run(
                ["Rscript", str(PLOT_R), str(d / "cov.bedgraph"), str(d / "cyto.txt"),
                 str(d / "out.pdf"), "NC_017628.1", "2193827", "2244594"],
                capture_output=True, text=True)
            for line in result.stderr.splitlines():
                if "coverage bins in" in line:
                    return int(line.split(":")[1].strip().split()[0])
            self.fail(f"no bin count reported. stderr:\n{result.stderr[-800:]}")

    def test_a_bin_straddling_the_island_start_is_kept(self):
        # Wholly outside [2193827, ...] by the old containment rule, so previously lost.
        self.assertEqual(self.bins_reported("NC_017628.1\t2193800\t2194000\t42\n"), 1)

    def test_a_bin_straddling_the_island_end_is_kept(self):
        self.assertEqual(self.bins_reported("NC_017628.1\t2244500\t2244700\t7\n"), 1)

    def test_a_bin_entirely_outside_the_island_is_still_dropped(self):
        self.assertEqual(self.bins_reported("NC_017628.1\t2100000\t2100500\t9\n"), 0)

    def test_a_bin_on_another_contig_is_dropped(self):
        self.assertEqual(self.bins_reported("chr1\t2193900\t2194000\t42\n"), 0)

    def test_interior_bins_are_unaffected(self):
        self.assertEqual(self.bins_reported(
            "NC_017628.1\t2200000\t2200500\t5\nNC_017628.1\t2210000\t2210500\t6\n"), 2)


if __name__ == "__main__":
    unittest.main()
