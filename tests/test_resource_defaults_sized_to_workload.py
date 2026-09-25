"""F17 (Ludmil, GATE 5): resource defaults sized from the measured workload, not guessed.

conf/base.config's per-process cpu/memory floors were, before commit e6d91f7, one to
three orders of magnitude above the peak RSS actually observed in the v0.0.1
Mutographs CRC trace (pksProfiler_analysis/v0.0.1/mutographs_crc/trace.txt):

    process            measured max      was asking
    extractReads           0.04 GB          64 GB
    filterReads             1.3 GB         128 GB
    mapReads               11.6 GB          64 GB
    pksProfilerAlign        1.7 GB         128 GB

On a shared cluster that costs throughput, not money: a 256 GB node can host two
filterReads tasks asking for 128 GB each instead of dozens asking for what they need.
No test asserted this stayed fixed -- a config edit could silently put any of these
back to its old, grossly oversized request and nothing here would notice.

This does not pin the exact current numbers (a legitimate future re-measurement should
not have to touch a test file to raise a floor a little), only that they stay far
below the specific pre-fix values this finding was about. Real conda-environment
execution coverage for these processes' actual behaviour lives in
tests/check_e2e_execution.sh; this is the static guard for the sizing itself.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_CONFIG = (ROOT / "conf/base.config").read_text()

# (label, the exact pre-fix GB request this finding cites, a generous ceiling well
# above the measured max but comfortably below the old value)
PRE_FIX_MEMORY_GB = {
    "extract_reads": (64, 32),
    "filter_reads": (128, 64),
    "map_reads": (64, 32),
    "pks_align": (128, 64),
}


def withlabel_block(label):
    match = re.search(
        r"withLabel:" + re.escape(label) + r"\s*\{(.*?)\n\s*\}",
        BASE_CONFIG, re.S)
    if not match:
        raise AssertionError(f"no withLabel:{label} block found in conf/base.config")
    return match.group(1)


def memory_gb(block):
    match = re.search(r"memory\s*=\s*\{?\s*(\d+(?:\.\d+)?)\s*\.GB", block)
    if not match:
        raise AssertionError(f"no memory = N.GB found in block: {block!r}")
    return float(match.group(1))


class ResourceFloorsStayFarBelowTheirPreFixValues(unittest.TestCase):
    def test_every_measured_label_is_below_its_pre_fix_request(self):
        for label, (pre_fix, ceiling) in PRE_FIX_MEMORY_GB.items():
            block = withlabel_block(label)
            current = memory_gb(block)
            self.assertLess(
                current, ceiling,
                f"withLabel:{label} requests {current} GB, no longer far below its "
                f"pre-fix {pre_fix} GB (F17: these were one to three orders of "
                "magnitude over the measured peak)")

    def test_map_reads_keeps_the_largest_floor_for_the_minimap2_index(self):
        # conf/base.config's own stated reason map_reads stays higher than the others:
        # minimap2 loads an 8.6 GB index. A future edit that quietly drops this below
        # the other per-sample labels would under-provision the one process that
        # actually needs the headroom.
        map_reads = memory_gb(withlabel_block("map_reads"))
        for other in ("extract_reads", "filter_reads", "pks_align"):
            self.assertGreaterEqual(map_reads, memory_gb(withlabel_block(other)))


if __name__ == "__main__":
    unittest.main()
