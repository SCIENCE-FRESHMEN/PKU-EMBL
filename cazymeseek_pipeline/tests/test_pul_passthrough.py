import tempfile
import unittest
from pathlib import Path

from cazymeseek_pipeline.standardize import standardize


class PulPassthroughTest(unittest.TestCase):
    def test_pul_homology_and_vote_are_exported_separately(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "abundance.tsv").write_text("protein_sequence_id\tTPM\np1\t10\n", encoding="utf-8")
            (root / "hmmer.tsv").write_text("Gene ID\tHMM Profile\tAli From\tAli To\np1\tGH1\t1\t50\n", encoding="utf-8")
            (root / "cgc.tsv").write_text("CGC ID\tGene ID\tContig ID\tAnnotation\nc1\tp1\tcontig\tGH1\n", encoding="utf-8")
            (root / "substrate.tsv").write_text("CGC ID\tSubstrate of the hit PUL\tSubstrate predicted by majority voting of CAZymes in CGC\nc1\txylan\tcellulose\n", encoding="utf-8")
            output = root / "out.csv"; standardize("s", root / "abundance.tsv", output, hmmer=root / "hmmer.tsv", cgc_standard=root / "cgc.tsv", substrate=root / "substrate.tsv")
            header, row = output.read_text(encoding="utf-8").splitlines(); columns = header.split(","); values = row.split(",")
            self.assertEqual(values[columns.index("CGC_substrate_PUL")], "xylan")
            self.assertEqual(values[columns.index("CGC_substrate_vote")], "cellulose")
