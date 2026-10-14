"""The README's Software table states exactly the versions pinned in conda_envs/*.yml."""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text()
PIN = re.compile(r"^\s*-\s+([A-Za-z0-9_.\-]+)\s*={1,2}\s*([^\s#]+)", re.M)


class ReadmeSoftwareTable(unittest.TestCase):
    def test_every_env_and_pin_is_in_the_table(self):
        section = README.split("## Software", 1)[1].split("\n## ", 1)[0]
        rows = dict(re.findall(r"^\| `([^`]+\.yml)` \| (.*) \|$", section, re.M))
        envs = sorted(p.name for p in (ROOT / "conda_envs").glob("*.yml"))
        self.assertEqual(sorted(rows), envs)
        for env in envs:
            pins = {f"{name} {version}" for name, version in PIN.findall((ROOT / "conda_envs" / env).read_text())}
            with self.subTest(env=env):
                self.assertEqual(set(rows[env].split(", ")), pins)


if __name__ == "__main__":
    unittest.main()
