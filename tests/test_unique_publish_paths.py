"""No two processes publish the same path.

plotPKSTaxa (--pks_taxa) and plotCommunityTaxa (--pks_community_taxa) are one process,
plotBrackenTaxa, included under two aliases, and both used to publish
by_sample/<s>/figures/taxa_barplots.pdf, so a run with both kept whichever finished last
(fixed in 6a03115: `<scope>_taxa_barplots.pdf`). Nextflow does not warn when two tasks
publish onto the same path; the second copy wins. tests/test_output_layout.py pins that
one case. This test checks every publishDir in Modules/ against every other, and on its
first run found the same overwrite one stage earlier: Bracken (--pks_taxa) and
sampleBracken (--pks_community_taxa) both published taxonomy/bracken.{G,S}.report.txt
and .krakenreport.txt (rc1.1: sampleBracken's are now community_bracken.*).

For each process it works out the published path templates: each publishDir's directory,
joined with each declared output name that the directive's `pattern` admits, renamed by
its `saveAs`. Interpolations stay as placeholders (`<sampleID>`, `<label>`), so two
templates collide when they are equal, or when one is a glob that matches the other. A
process included under more than one name collides with itself unless its published
paths carry a `val` input, and every call site passes that input a different literal.

saveAs forms understood: `{ fn -> fn - "X" }`, `{ "X" }` and no saveAs. Any other form
must be added to SAVE_AS below with the rename it performs, or this test fails, so a new
rename is never silently skipped.
"""
import fnmatch
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_FILES = sorted(p for p in (ROOT / "Modules").glob("*.nf") if ".bak" not in p.name)
MAIN = (ROOT / "main.nf").read_text()

PROCESS = re.compile(r"^process\s+(\w+)\s*\{", re.M)
SECTION = re.compile(r"^\s*(input|output|script|shell|exec|stub|when)\s*:", re.M)


def template(text):
    """`${sampleID}.x` and `$sampleID.x` -> `<sampleID>.x`; `${params.sample_dir}` -> `<params.sample_dir>`."""
    text = re.sub(r"\$\{\s*([^}]+?)\s*\}", r"<\1>", text)
    return re.sub(r"\$(\w+)", r"<\1>", text)


def concrete(path):
    """A template as a filename: every placeholder becomes one fixed token."""
    return re.sub(r"<[^>]+>", "VAR", path)


# saveAs closures other than the two common forms, with the rename each performs.
def _align_save_as(name):
    # pksProfiler_align: counts and coverage at the sample root, the rest under alignment/.
    stripped = name.replace("<sampleID>.", "", 1)
    if name.endswith(".counts.txt") or name.endswith(".coverage.txt"):
        return stripped
    return "alignment/" + stripped


SAVE_AS = {
    # sampleBracken (--pks_community_taxa): distinct from Bracken's (--pks_taxa) names
    'fn -> "community_" + (fn - "${sampleID}.")':
        lambda name: "community_" + name.replace("<sampleID>.", "", 1),
    # a directory output published flat: its basename
    "fn -> fn.tokenize('/').last()": lambda name: name.rstrip("/").rsplit("/", 1)[-1],
    "fn -> (fn.endsWith('.counts.txt') || fn.endsWith('.coverage.txt')) "
    "? fn - \"${sampleID}.\" : \"alignment/\" + (fn - \"${sampleID}.\")": _align_save_as,
}


def balanced(text, start):
    """The text between the bracket at `start` and its match (Groovy, minus strings)."""
    pairs = {"(": ")", "{": "}"}
    stack, i, quote = [], start, None
    while i < len(text):
        c = text[i]
        if quote:
            if c == "\\":
                i += 2
                continue
            if c == quote:
                quote = None
            elif quote == '"' and text.startswith("${", i):
                # a GString interpolation: skip to its closing brace
                i = text.index("}", i)
        elif c in "\"'":
            quote = c
        elif c in pairs:
            stack.append(pairs[c])
        elif stack and c == stack[-1]:
            stack.pop()
            if not stack:
                return text[start + 1:i]
        i += 1
    raise ValueError(f"unbalanced bracket at {start}")


def strip_comment(line):
    """A Groovy line without its `//` comment (a `//` inside a string is kept)."""
    quote = None
    for i, c in enumerate(line):
        if quote:
            if c == quote:
                quote = None
        elif c in "\"'":
            quote = c
        elif line.startswith("//", i):
            return line[:i]
    return line


def publish_directives(header):
    """[(dir, pattern or None, saveAs closure body or None)] for one process's directives."""
    out = []
    for match in re.finditer(r"\bpublishDir\b\s*(\(?)", header):
        if match.group(1):
            args = balanced(header, match.end() - 1)
        else:
            # Unparenthesised: the directive runs to the end of its logical line, which
            # continues while a closure is open.
            start, depth, i = match.end(), 0, match.end()
            while i < len(header):
                if header[i] == "{":
                    body = balanced(header, i)
                    i += len(body) + 2
                    continue
                if header[i] == "\n" and not header[start:i].rstrip().endswith(","):
                    break
                i += 1
            args = header[start:i]
        if args.lstrip().startswith("{"):
            first = args.index("{")
            directory = balanced(args, first).strip().strip('"')
            rest = args[first + len(balanced(args, first)) + 2:]
        else:
            quoted = re.match(r'\s*"([^"]*)"', args)
            directory, rest = quoted.group(1), args[quoted.end():]
        pattern = re.search(r'pattern:\s*"([^"]*)"', rest)
        save_as = None
        at = rest.find("saveAs:")
        if at >= 0:
            brace = rest.index("{", at)
            save_as = " ".join(balanced(rest, brace).split())
        out.append((template(directory), pattern.group(1) if pattern else None, save_as))
    return out


def output_names(body):
    """Templates of the declared output files, and the bare identifiers among them."""
    sections = list(SECTION.finditer(body))
    names, bare = [], []
    for i, sec in enumerate(sections):
        if sec.group(1) != "output":
            continue
        end = sections[i + 1].start() if i + 1 < len(sections) else len(body)
        block = "\n".join(strip_comment(l) for l in body[sec.end():end].splitlines())
        for m in re.finditer(r"\bpath\s*\(?\s*(\"([^\"]*)\"|'([^']*)'|(\w+))", block):
            if m.group(4):
                bare.append(m.group(4))
            else:
                names.append(template(m.group(2) if m.group(2) is not None else m.group(3)))
    return names, bare


def input_names(body):
    sections = list(SECTION.finditer(body))
    for i, sec in enumerate(sections):
        if sec.group(1) == "input":
            end = sections[i + 1].start() if i + 1 < len(sections) else len(body)
            block = body[sec.end():end]
            return (set(re.findall(r"\bpath\s*\(?\s*(\w+)", block)),
                    re.findall(r"\bval\s*\(?\s*(\w+)", block))
    return set(), []


def rename(save_as, name):
    if save_as is None:
        return name
    constant = re.fullmatch(r'"([^"]*)"', save_as)
    if constant:
        return template(constant.group(1))
    minus = re.fullmatch(r'fn\s*->\s*fn\s*-\s*"([^"]*)"', save_as)
    if minus:
        return name.replace(template(minus.group(1)), "", 1)
    if save_as in SAVE_AS:
        return SAVE_AS[save_as](name)
    raise AssertionError(f"saveAs form not understood, add it to SAVE_AS: {{ {save_as} }}")


UNKNOWN_SAVE_AS = []


def processes():
    """{process name: {"file", "published": [templates], "vals": [...], "bare": [...]}}."""
    out = {}
    for path in MODULE_FILES:
        text = path.read_text()
        starts = list(PROCESS.finditer(text))
        for index, match in enumerate(starts):
            # Up to the next process: the script block is shell, and its quotes and
            # braces are not Groovy's, so the body is not bracket-matched.
            end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
            body = text[match.end():end]
            first = SECTION.search(body)
            header = "\n".join(strip_comment(l) for l in
                               (body[:first.start()] if first else body).splitlines())
            names, bare = output_names(body)
            inputs, vals = input_names(body)
            published = []
            for directory, pattern, save_as in publish_directives(header):
                for name in names:
                    if pattern and not fnmatch.fnmatch(concrete(name), pattern):
                        continue
                    try:
                        published.append(f"{directory}/{rename(save_as, name)}")
                    except AssertionError as error:
                        UNKNOWN_SAVE_AS.append(f"{path.name}:{match.group(1)}: {error}")
            out[match.group(1)] = {"file": path.name, "published": published, "vals": vals,
                                   "bare": [b for b in bare if b not in inputs],
                                   "has_publish": "publishDir" in header}
    return out


def collide(a, b):
    if a == b:
        return True
    ca, cb = concrete(a), concrete(b)
    return ("*" in a and fnmatch.fnmatch(cb, ca)) or ("*" in b and fnmatch.fnmatch(ca, cb))


def include_names():
    """{declared process: [names it is included under]} across main.nf and Modules/."""
    names = {}
    for text in [MAIN] + [p.read_text() for p in MODULE_FILES]:
        for block in re.findall(r"include\s*\{([^}]*)\}", text):
            for item in block.split(";"):
                parts = item.split()
                if not parts:
                    continue
                declared = parts[0]
                alias = parts[2] if len(parts) == 3 and parts[1] == "as" else declared
                names.setdefault(declared, [])
                if alias not in names[declared]:
                    names[declared].append(alias)
    return names


PROCESSES = processes()


class EveryPublishedPathHasOneWriter(unittest.TestCase):
    def test_every_save_as_is_understood(self):
        self.assertEqual(UNKNOWN_SAVE_AS, [])

    def test_the_parser_sees_every_publishing_process(self):
        # A process whose publishDir the parser cannot resolve to any path would be
        # silently exempt.
        empty = [n for n, p in PROCESSES.items() if p["has_publish"] and not p["published"]]
        self.assertEqual(empty, [])
        self.assertGreater(len(PROCESSES), 40)
        self.assertGreater(sum(len(p["published"]) for p in PROCESSES.values()), 60)

    def test_every_output_name_is_a_literal_or_an_input(self):
        # A bare identifier that is not a staged input is a name the test cannot see.
        bare = {n: p["bare"] for n, p in PROCESSES.items() if p["has_publish"] and p["bare"]}
        self.assertEqual(bare, {})

    def test_no_two_processes_publish_the_same_path(self):
        clashes = []
        items = [(name, path) for name, p in PROCESSES.items() for path in set(p["published"])]
        for i, (name_a, path_a) in enumerate(items):
            for name_b, path_b in items[i + 1:]:
                if name_a != name_b and collide(path_a, path_b):
                    clashes.append(f"{name_a} and {name_b}: {path_a} / {path_b}")
        self.assertEqual(clashes, [])

    def test_a_process_included_twice_publishes_under_a_name_its_call_site_sets(self):
        for declared, aliases in include_names().items():
            if len(aliases) < 2 or declared not in PROCESSES:
                continue
            proc = PROCESSES[declared]
            if not proc["published"]:
                continue
            with self.subTest(process=declared, aliases=aliases):
                keys = [v for v in proc["vals"] if v != "sampleID"]
                distinguishing = [v for v in keys
                                  if all(f"<{v}>" in path for path in proc["published"])]
                self.assertTrue(distinguishing,
                                f"{declared} is used as {aliases} but every alias publishes "
                                f"the same paths: {sorted(set(proc['published']))}")
                literals = {}
                for alias in aliases:
                    calls = re.findall(rf"\b{alias}\s*\((.*)\)\s*$", MAIN, re.M)
                    self.assertEqual(len(calls), 1, f"{alias}: expected one call site")
                    literal = re.findall(r"'([^']*)'|\"([^\"]*)\"", calls[0])
                    self.assertTrue(literal, f"{alias}: the call passes no literal to tell it apart")
                    literals[alias] = "".join(literal[-1])
                self.assertEqual(len(set(literals.values())), len(literals), literals)

    def test_the_two_bracken_runs_publish_distinct_tables(self):
        # Found by this test: the same overwrite, one stage before the barplots.
        island = set(PROCESSES["Bracken"]["published"])
        community = set(PROCESSES["sampleBracken"]["published"])
        self.assertIn("<params.sample_dir>/<sampleID>/taxonomy/bracken.S.report.txt", island)
        self.assertIn("<params.sample_dir>/<sampleID>/taxonomy/community_bracken.S.report.txt", community)
        self.assertEqual(island & community, set())

    def test_the_two_taxa_barplots_are_distinct(self):
        # The case that motivated this file.
        published = PROCESSES["plotBrackenTaxa"]["published"]
        self.assertEqual(published, ["<params.sample_dir>/<sampleID>/figures/<scope>_taxa_barplots.pdf"])
        self.assertIn("plotPKSTaxa(BRACKEN_GS_REPORTS, 'pks_island')", MAIN)
        self.assertRegex(MAIN, r"plotCommunityTaxa\(.*, 'community'\)")


class TheCheckCatchesACollision(unittest.TestCase):
    """The collision rules themselves, on the shapes that occur."""

    def test_equal_templates_collide(self):
        self.assertTrue(collide("<d>/<sampleID>/figures/taxa_barplots.pdf",
                                "<d>/<sampleID>/figures/taxa_barplots.pdf"))

    def test_a_glob_collides_with_what_it_matches(self):
        self.assertTrue(collide("<d>/<sampleID>/figures/*.pdf", "<d>/<sampleID>/figures/circos.pdf"))
        self.assertFalse(collide("<d>/<sampleID>/figures/*.svg", "<d>/<sampleID>/figures/circos.pdf"))

    def test_different_names_or_dirs_do_not(self):
        self.assertFalse(collide("<d>/<sampleID>/figures/pks_island_taxa_barplots.pdf",
                                 "<d>/<sampleID>/figures/community_taxa_barplots.pdf"))
        self.assertFalse(collide("<d>/<sampleID>/hmm/x.tsv", "<d>/<sampleID>/taxonomy/x.tsv"))

    def test_rename_forms(self):
        self.assertEqual(rename('fn -> fn - "${sampleID}."', "<sampleID>.counts.txt"), "counts.txt")
        self.assertEqual(rename('"${scope}_taxa_barplots.pdf"', "<sampleID>.x.pdf"),
                         "<scope>_taxa_barplots.pdf")
        self.assertEqual(rename(None, "<sampleID>.x"), "<sampleID>.x")
        with self.assertRaises(AssertionError):
            rename("fn -> fn.toUpperCase()", "x")

    def test_the_old_shared_name_would_have_been_caught(self):
        old = rename('"taxa_barplots.pdf"', "<sampleID>.pks_island_taxa_barplots.pdf")
        self.assertNotIn("<scope>", old)


if __name__ == "__main__":
    unittest.main()
