import unittest

from cazymeseek_pipeline.annotation import ec_values


class EcMappingTest(unittest.TestCase):
    def test_only_dbcan_sub_domain_receives_curated_ec(self):
        overview = {"EC#": "3.2.1.21"}
        self.assertEqual(ec_values(overview, {"method": "dbCAN_sub", "domain_ec": ""}), "3.2.1.21")
        self.assertEqual(ec_values(overview, {"method": "dbCAN", "domain_ec": ""}), "")

    def test_missing_mapping_stays_blank(self):
        self.assertEqual(ec_values({}, {"method": "dbCAN_sub", "domain_ec": ""}), "")
