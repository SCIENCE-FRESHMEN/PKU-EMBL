import tempfile
import unittest
from pathlib import Path

from cazymeseek_pipeline.deduplicate import read_clusters, write_representative_fasta


class DeduplicationTest(unittest.TestCase):
    def test_cluster_mapping_and_representative_cds_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); clusters = root / "clusters.tsv"; fasta = root / "genes.ffn"; output = root / "nr.ffn"
            clusters.write_text("p1\tp1\np1\tp2\n", encoding="utf-8")
            fasta.write_text(">p1\nATGC\n>p2\nGGCC\n", encoding="utf-8")
            self.assertEqual(read_clusters(clusters), {"p1": "p1", "p2": "p1"})
            write_representative_fasta(fasta, clusters, output)
            self.assertEqual(output.read_text(encoding="utf-8"), ">p1\nATGC\n")
