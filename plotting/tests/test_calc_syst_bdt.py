"""Regression tests for the two-sided systematic aggregation; run in the PPG12 ROOT environment.

Convention under test (decision 2026-09-08, matches the paper systematics):
the "down" variant fills the low column and the "up" variant fills the high
column with |variant - nominal|, whatever direction the spectrum moved. A
missing direction is mirrored from the other one. The alternative signed
per-bin envelope is preserved in reports/two_sided_fix_2026-09-08/ and is
not what the pipeline uses.
"""
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import calc_syst_bdt as syst


def histogram(name, values):
    h = ROOT.TH1D(name, "", len(values), 0., float(len(values)))
    h.SetDirectory(0)
    for i, value in enumerate(values, 1):
        h.SetBinContent(i, value)
    return h


class TwoSidedEnvelopeTest(unittest.TestCase):
    def aggregate(self, nominal, down, up, skip_missing=False):
        h_nom = histogram("test_nominal", nominal)
        inputs = {}
        for role, values in [("down", down), ("up", up)]:
            if values is not None:
                inputs["bdt_test_" + role] = histogram("test_" + role, values)

        def load(name, *_args):
            if name not in inputs:
                raise FileNotFoundError(name)
            return inputs[name]

        vmap = {"npb_cut": {"down": ["test_down"], "up": ["test_up"],
                            "one_sided": [], "max": []}}
        with patch.object(syst, "load_spectrum", side_effect=load):
            return syst.aggregate_type("npb_cut", vmap, h_nom, "", "", skip_missing)

    def assert_contents(self, h, expected):
        self.assertEqual(h.GetNbinsX(), len(expected))
        for i, value in enumerate(expected, 1):
            with self.subTest(bin=i):
                self.assertAlmostEqual(h.GetBinContent(i), value, places=12)

    def test_role_based_columns_and_relative_errors(self):
        # Down variant -> low column, up variant -> high column, as |var - nom|,
        # independent of the direction the spectrum actually moved.
        result = self.aggregate(
            [100.] * 6,
            [105., 94., 103., 92., 100., 100.],
            [90., 108., 107., 97., 100., 104.],
        )
        for h, expected in zip(result, [
            [5., 6., 3., 8., 0., 0.],
            [10., 8., 7., 3., 0., 4.],
            [.05, .06, .03, .08, 0., 0.],
            [.10, .08, .07, .03, 0., .04],
        ]):
            self.assert_contents(h, expected)

    def test_swapping_parameter_labels_swaps_columns(self):
        a = self.aggregate([100., 200.], [105., 180.], [90., 206.])
        b = self.aggregate([100., 200.], [90., 206.], [105., 180.])
        # (low, high, rel_low, rel_high) of a == (high, low, rel_high, rel_low) of b
        for ha, hb in zip(a, (b[1], b[0], b[3], b[2])):
            self.assert_contents(ha, [hb.GetBinContent(i) for i in (1, 2)])

    def test_zero_and_negative_nominal(self):
        result = self.aggregate([0., -100.], [-2., -110.], [3., -95.])
        for h, expected in zip(result, [[2., 10.], [3., 5.], [0., .10], [0., .05]]):
            self.assert_contents(h, expected)

    def test_missing_direction_keeps_symmetric_fallback(self):
        for down, up in [(None, [110.]), ([90.], None)]:
            with self.subTest(down=down, up=up):
                result = self.aggregate([100.], down, up, skip_missing=True)
                for h, expected in zip(result, [[10.], [10.], [.10], [.10]]):
                    self.assert_contents(h, expected)

    def test_missing_required_direction_still_raises_in_aggregator(self):
        with self.assertRaises(FileNotFoundError):
            self.aggregate([100.], [90.], None)

    def test_role_based_columns_reach_total_root_file(self):
        # Two independent sources with down-variant shifts 3 and 3 and
        # up-variant shifts 4 and 4; the first source moves opposite to its labels.
        first = self.aggregate([100.], [103.], [96.])
        second = self.aggregate([100.], [97.], [104.])
        total = syst.quadrature_sum([first, second])
        low, high = math.sqrt(3.**2 + 3.**2), math.sqrt(4.**2 + 4.**2)
        with tempfile.TemporaryDirectory(prefix="ppg12-syst-test-") as directory:
            syst.write_syst_sum(directory, *total)
            f = ROOT.TFile.Open(str(Path(directory) / "syst_sum.root"))
            try:
                for name, value in [
                    ("h_sum_low", low), ("h_sum_high", high),
                    ("h_sum_rel_low", low / 100.), ("h_sum_rel_high", high / 100.),
                ]:
                    self.assert_contents(f.Get(name), [value])
            finally:
                f.Close()


if __name__ == "__main__":
    unittest.main()
