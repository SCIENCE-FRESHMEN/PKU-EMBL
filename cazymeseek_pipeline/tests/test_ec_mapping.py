import unittest

from cazymeseek_pipeline.annotation import map_ec


class EcMappingTest(unittest.TestCase):
    def test_subfamily_mapping_precedes_family_mapping(self):
        candidate = {"cazy_family": "GH43", "cazy_subfamily": "GH43_4"}
        self.assertEqual(map_ec(candidate, {"GH43_4": "3.2.1.55", "GH43": "3.2.1.99"}), "3.2.1.55")

    def test_family_is_fallback_and_missing_is_blank(self):
        candidate = {"cazy_family": "GH43", "cazy_subfamily": "GH43_4"}
        self.assertEqual(map_ec(candidate, {"GH43": "3.2.1.99"}), "3.2.1.99")
        self.assertEqual(map_ec(candidate, {}), "")
