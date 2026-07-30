import unittest

from cazymeseek_pipeline.annotation import preferred_family


class ConflictPriorityTest(unittest.TestCase):
    def test_dbcan_hmm_overrides_sub_and_diamond_conflict(self):
        family, method, calls = preferred_family({"dbCAN": "GH1(5-90)", "dbCAN_sub": "GH43_e149", "DIAMOND": "GH43_4"})
        self.assertEqual((family, method), ("GH1", "dbCAN"))
        self.assertEqual(calls["DIAMOND"], ["GH43_4"])

    def test_dbcan_sub_overrides_diamond_when_hmm_is_absent(self):
        family, method, _ = preferred_family({"dbCAN": "-", "dbCAN_sub": "GH43_e149", "DIAMOND": "GH43_4"})
        self.assertEqual((family, method), ("GH43", "dbCAN_sub"))
