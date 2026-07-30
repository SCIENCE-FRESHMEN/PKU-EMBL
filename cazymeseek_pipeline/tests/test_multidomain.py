import unittest

from cazymeseek_pipeline.annotation import resolve_conflicts


class MultiDomainTest(unittest.TestCase):
    def test_non_overlapping_top_rank_domains_are_all_retained(self):
        calls = [
            {"gene_id": "p1", "cazy_family": "GH1", "cazy_subfamily": "GH1", "domain_start": 1, "domain_end": 80, "annotation_source": "hmmer", "source_rank": 0},
            {"gene_id": "p1", "cazy_family": "GH3", "cazy_subfamily": "GH3", "domain_start": 120, "domain_end": 220, "annotation_source": "hmmer", "source_rank": 0},
            {"gene_id": "p1", "cazy_family": "GH43", "cazy_subfamily": "GH43_4", "domain_start": 130, "domain_end": 215, "annotation_source": "diamond", "source_rank": 2},
        ]
        kept = resolve_conflicts(calls, .80)
        self.assertEqual([(item["cazy_family"], item["domain_start"]) for item in kept], [("GH1", 1), ("GH3", 120)])
