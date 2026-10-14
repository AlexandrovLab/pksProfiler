"""F17 (Ludmil, 19 and 23 Sep reports): the Biowulf profile must use NIH's documented
scheduler settings, queueStatInterval 5 min and submitRateLimit 6/1min. rc1 reported
F17 fixed but left these at 2 min and 1 sec."""
import re
import unittest
from pathlib import Path

CONFIG = (Path(__file__).resolve().parents[1] / "conf/biowulf.config").read_text()


def setting(name):
    match = re.search(rf"^\s*{name}\s*=\s*'([^']+)'", CONFIG, re.M)
    return match.group(1) if match else None


class BiowulfScheduler(unittest.TestCase):
    def test_queue_stat_interval(self):
        self.assertEqual(setting("queueStatInterval"), "5 min")

    def test_submit_rate_limit(self):
        self.assertEqual(setting("submitRateLimit"), "6/1min")

    def test_poll_interval_is_nih_default(self):
        self.assertEqual(setting("pollInterval"), "2 min")


if __name__ == "__main__":
    unittest.main()
