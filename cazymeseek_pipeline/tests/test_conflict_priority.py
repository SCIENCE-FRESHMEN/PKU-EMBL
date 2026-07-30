import csv
import tempfile
import unittest
from pathlib import Path

from cazymeseek_pipeline.annotation import resolve_conflicts
from cazymeseek_pipeline.standardize import standardize


class ConflictPriorityTest(unittest.TestCase):
    def test_hmmer_wins_over_hotpep_and_diamond_at_same_domain(self):
        calls = [
            {"gene_id": "p1", "cazy_family": "GH1", "cazy_subfamily": "GH1", "domain_start": 10, "domain_end": 100, "annotation_source": "hmmer", "source_rank": 0},
            {"gene_id": "p1", "cazy_family": "GH43", "cazy_subfamily": "GH43_4", "domain_start": 12, "domain_end": 98, "annotation_source": "hotpep", "source_rank": 1},
            {"gene_id": "p1", "cazy_family": "GH43", "cazy_subfamily": "GH43_4", "domain_start": 10, "domain_end": 100, "annotation_source": "diamond", "source_rank": 2},
        ]
        kept = resolve_conflicts(calls, .80)
        self.assertEqual([(item["annotation_source"], item["source_rank"] ) for item in kept], [("hmmer", 0)])

    def test_standardized_output_contains_source_rank(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "hmmer.tsv").write_text("Gene ID\tHMM Profile\tAli From\tAli To\np1\tGH1\t1\t100\n", encoding="utf-8")
            (root / "hotpep.tsv").write_text("Gene ID\tHMM Profile\tAli From\tAli To\np1\tGH43_4\t2\t99\n", encoding="utf-8")
            (root / "abundance.tsv").write_text("protein_sequence_id\tTPM\np1\t10\n", encoding="utf-8")
            output = root / "out.csv"; standardize("s", root / "abundance.tsv", output, hmmer=root / "hmmer.tsv", hotpep=root / "hotpep.tsv")
            with output.open(newline="", encoding="utf-8") as handle: row = next(csv.DictReader(handle))
            self.assertEqual((row["annotation_source"], row["source_rank"]), ("hmmer", "0"))
