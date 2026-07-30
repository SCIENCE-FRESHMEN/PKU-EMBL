import tempfile
import unittest
from pathlib import Path

from cazymeseek_pipeline.standardize import standardize


class MultiDomainTest(unittest.TestCase):
    def test_all_hmm_domains_become_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "overview.txt").write_text("Gene ID\tEC#\tdbCAN\tdbCAN_sub\tDIAMOND\np1\t\tGH1;GH3\t-\t-\n", encoding="utf-8")
            (root / "hmmer.out").write_text("Gene ID\tHMM Profile\tAli From\tAli To\np1\tGH1\t1\t50\np1\tGH3\t60\t120\n", encoding="utf-8")
            (root / "abundance.tsv").write_text("protein_sequence_id\tTPM\np1\t10\n", encoding="utf-8")
            output = root / "out.csv"; standardize("s", root / "overview.txt", root / "abundance.tsv", output, root / "hmmer.out")
            rows = output.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(rows), 3)
            self.assertIn(",GH1,GH1,", rows[1]); self.assertIn(",GH3,GH3,", rows[2])
