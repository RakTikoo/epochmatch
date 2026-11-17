"""Tests for epochmatch.

Run with:  python3 run_tests.py
"""

import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from epochmatch.cli import main as cli_main
from epochmatch.core import match_tables, propagate_positions


def make_table(ids, ra, dec, pmra=None, pmdec=None):
    t = Table()
    t["id"] = ids
    t["ra"] = ra
    t["dec"] = dec
    if pmra is not None:
        t["pmra"] = pmra
    if pmdec is not None:
        t["pmdec"] = pmdec
    return t


class TestEpochPropagationRecoversMatches(unittest.TestCase):
    """(1) Stars that moved far between epochs: a naive position match
    fails, epochmatch recovers the true pairs."""

    def setUp(self):
        # Epoch A = 2000.0.  Motions ~600 mas/yr over 15.5 yr move the
        # stars ~9 arcsec -- far outside a 2 arcsec match radius.
        self.a = make_table(
            ["a1", "a2", "a3"],
            [10.0, 10.05, 10.10],
            [30.0, 30.05, 30.10],
            pmra=[500.0, -450.0, 100.0],
            pmdec=[300.0, 200.0, -550.0],
        )
        prop = propagate_positions(self.a["ra"], self.a["dec"],
                                   self.a["pmra"], self.a["pmdec"],
                                   2000.0, 2015.5)
        self.b = make_table(
            ["b1", "b2", "b3"],
            list(prop.ra.deg), list(prop.dec.deg),
            pmra=[500.0, -450.0, 100.0],
            pmdec=[300.0, 200.0, -550.0],
        )

    def test_naive_static_match_would_fail(self):
        ca = SkyCoord(ra=self.a["ra"] * u.deg, dec=self.a["dec"] * u.deg)
        cb = SkyCoord(ra=self.b["ra"] * u.deg, dec=self.b["dec"] * u.deg)
        _, sep, _ = ca.match_to_catalog_sky(cb)
        self.assertTrue(np.all(sep.arcsec > 2.0),
                        "fixture broken: stars did not move enough")

    def test_epochmatch_recovers_true_pairs(self):
        result = match_tables(self.a, self.b, epoch_a=2000.0,
                              epoch_b=2015.5, to_epoch=2016.0,
                              radius_arcsec=2.0, id_a="id", id_b="id")
        pairs = {(r["id_a"], r["id_b"]) for r in result.matched}
        self.assertEqual(pairs, {("a1", "b1"), ("a2", "b2"), ("a3", "b3")})
        self.assertTrue(np.all(result.matched["sep_arcsec"] < 0.01))
        self.assertEqual(len(result.unmatched_a), 0)
        self.assertEqual(len(result.unmatched_b), 0)
        self.assertEqual(len(result.ambiguous), 0)


class TestZeroMotionEqualsStaticMatch(unittest.TestCase):
    """(2) With zero proper motion everywhere, epochmatch must agree
    with a plain static SkyCoord match."""

    def test_zero_motion(self):
        a = make_table(["a1", "a2", "a3"], [10.0, 10.05, 10.10],
                       [30.0, 30.05, 30.10], pmra=[0.0] * 3,
                       pmdec=[0.0] * 3)
        # B rows offset by 0.5 arcsec in dec from A rows 1 and 3 only.
        b = make_table(["b1", "b3"], [10.0, 10.10],
                       [30.0 + 0.5 / 3600.0, 30.10 + 0.5 / 3600.0],
                       pmra=[0.0] * 2, pmdec=[0.0] * 2)
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2015.5,
                              to_epoch=2016.0, radius_arcsec=2.0,
                              id_a="id", id_b="id")
        ca = SkyCoord(ra=a["ra"] * u.deg, dec=a["dec"] * u.deg)
        cb = SkyCoord(ra=b["ra"] * u.deg, dec=b["dec"] * u.deg)
        idx, sep, _ = ca.match_to_catalog_sky(cb)
        expected = {(a["id"][i], b["id"][idx[i]])
                    for i in range(len(a)) if sep[i].arcsec <= 2.0}
        got = {(r["id_a"], r["id_b"]) for r in result.matched}
        self.assertEqual(got, expected)
        self.assertEqual(got, {("a1", "b1"), ("a3", "b3")})
        self.assertEqual([r["id_a"] for r in result.unmatched_a], ["a2"])


class TestAmbiguous(unittest.TestCase):
    """(3) Two B candidates inside the radius -> ambiguous, not matched."""

    def test_ambiguous_pair(self):
        a = make_table(["a1"], [10.0], [30.0], pmra=[0.0], pmdec=[0.0])
        b = make_table(["b1", "b2"],
                       [10.0, 10.0 + 1.0 / 3600.0],
                       [30.0 + 1.0 / 3600.0, 30.0],
                       pmra=[0.0, 0.0], pmdec=[0.0, 0.0])
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2000.0,
                              to_epoch=2000.0, radius_arcsec=2.0,
                              id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 0)
        self.assertEqual(len(result.ambiguous), 2)
        self.assertEqual({r["id_b"] for r in result.ambiguous},
                         {"b1", "b2"})
        self.assertEqual(len(result.unmatched_a), 0)
        self.assertEqual(len(result.unmatched_b), 0)


class TestMissingProperMotion(unittest.TestCase):
    """(4) A catalogue without PM columns warns and still matches
    statically."""

    def test_missing_pm_warns_and_works(self):
        a = make_table(["a1"], [10.0], [30.0], pmra=[0.0], pmdec=[0.0])
        b = make_table(["b1"], [10.0], [30.0])  # no pm columns at all
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2015.5,
                              to_epoch=2016.0, radius_arcsec=2.0,
                              id_a="id", id_b="id")
        self.assertTrue(any("no proper-motion" in w
                            for w in result.warnings))
        self.assertEqual(len(result.matched), 1)


class TestEdgeCases(unittest.TestCase):
    """Adversarial cases: empty input, NaN coordinates, duplicates,
    radius 0."""

    def test_empty_catalogue(self):
        a = make_table([], [], [], pmra=[], pmdec=[])
        b = make_table(["b1"], [10.0], [30.0], pmra=[0.0], pmdec=[0.0])
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2015.5,
                              to_epoch=2016.0, radius_arcsec=2.0,
                              id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 0)
        self.assertEqual(len(result.unmatched_b), 1)

    def test_nan_coordinates_go_to_unmatched_with_warning(self):
        a = make_table(["a1", "a2"], [10.0, np.nan], [30.0, 30.0],
                       pmra=[0.0, 0.0], pmdec=[0.0, 0.0])
        b = make_table(["b1"], [10.0], [30.0], pmra=[0.0], pmdec=[0.0])
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2015.5,
                              to_epoch=2016.0, radius_arcsec=2.0,
                              id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 1)
        self.assertEqual([r["id_a"] for r in result.unmatched_a], ["a2"])
        self.assertTrue(any("non-finite" in w for w in result.warnings))

    def test_duplicate_positions_are_ambiguous(self):
        a = make_table(["a1"], [10.0], [30.0], pmra=[0.0], pmdec=[0.0])
        b = make_table(["b1", "b2"], [10.0, 10.0], [30.0, 30.0],
                       pmra=[0.0, 0.0], pmdec=[0.0, 0.0])
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2000.0,
                              to_epoch=2000.0, radius_arcsec=2.0,
                              id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 0)
        self.assertEqual(len(result.ambiguous), 2)

    def test_radius_zero_matches_only_exact_coincidence(self):
        a = make_table(["a1", "a2"], [10.0, 10.01], [30.0, 30.0],
                       pmra=[0.0, 0.0], pmdec=[0.0, 0.0])
        b = make_table(["b1", "b2"], [10.0, 10.01],
                       [30.0, 30.0 + 0.1 / 3600.0],
                       pmra=[0.0, 0.0], pmdec=[0.0, 0.0])
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2000.0,
                              to_epoch=2000.0, radius_arcsec=0.0,
                              id_a="id", id_b="id")
        self.assertEqual({(r["id_a"], r["id_b"]) for r in result.matched},
                         {("a1", "b1")})


class TestCliErrors(unittest.TestCase):
    """User errors must be one clean line, never a traceback."""

    def _run_cli(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = cli_main(argv)
        return code, out.getvalue(), err.getvalue()

    def test_missing_file(self):
        code, _, err = self._run_cli(
            ["nope.csv", "alsono.fits", "--epoch-a", "2000",
             "--epoch-b", "2015.5"])
        self.assertEqual(code, 2)
        self.assertIn("file not found", err)
        self.assertNotIn("Traceback", err)

    def test_missing_column(self):
        with tempfile.TemporaryDirectory() as d:
            bad = os.path.join(d, "bad.csv")
            with open(bad, "w") as f:
                f.write("id,ra\nx,10.0\n")  # no dec column
            good = os.path.join(d, "good.csv")
            with open(good, "w") as f:
                f.write("ra,dec,pmra,pmdec\n10.0,30.0,0.0,0.0\n")
            code, _, err = self._run_cli(
                [bad, good, "--epoch-a", "2000", "--epoch-b", "2015.5"])
        self.assertEqual(code, 2)
        self.assertIn("missing Dec column", err)
        self.assertNotIn("Traceback", err)

    def test_cli_end_to_end_csv(self):
        with tempfile.TemporaryDirectory() as d:
            pa = os.path.join(d, "a.csv")
            pb = os.path.join(d, "b.csv")
            with open(pa, "w") as f:
                f.write("ra,dec,pmra,pmdec\n10.0,30.0,0.0,0.0\n")
            with open(pb, "w") as f:
                f.write("ra,dec,pmra,pmdec\n10.0,30.0,0.0,0.0\n")
            out_prefix = os.path.join(d, "res")
            code, out, _ = self._run_cli(
                [pa, pb, "--epoch-a", "2000", "--epoch-b", "2015.5",
                 "--out", out_prefix])
            self.assertEqual(code, 0)
            self.assertIn("Matched: 1", out)
            for kind in ("matched", "unmatched_a", "unmatched_b",
                         "ambiguous"):
                self.assertTrue(
                    os.path.isfile(f"{out_prefix}_{kind}.csv"), kind)


class TestAdversarialRegressions(unittest.TestCase):
    """Regressions from the independent adversarial review (REVIEW.md):
    non-finite/masked proper motions, non-finite radius/epochs, masked
    coordinates, radius-boundary and >180-degree radii, a B row claimed
    by both an ambiguous and a single-candidate A row, and Dec range."""

    def _one(self, **kw):
        return make_table(["x"], [10.0], [30.0],
                          pmra=[0.0], pmdec=[0.0], **kw)

    def test_nan_pm_in_a_goes_to_unmatched_with_warning(self):
        a = make_table(["a"], [10.0], [30.0], pmra=[np.nan], pmdec=[0.0])
        result = match_tables(a, self._one(), epoch_a=2000.0,
                              epoch_b=2015.5, to_epoch=2016.0,
                              radius_arcsec=2.0, id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 0)
        self.assertEqual([r["id_a"] for r in result.unmatched_a], ["a"])
        self.assertTrue(any("proper motion" in w and "non-finite" in w
                            for w in result.warnings))

    def test_nan_pm_in_b_goes_to_unmatched_with_warning(self):
        b = make_table(["b"], [10.0], [30.0], pmra=[0.0], pmdec=[np.nan])
        result = match_tables(self._one(), b, epoch_a=2000.0,
                              epoch_b=2015.5, to_epoch=2016.0,
                              radius_arcsec=2.0, id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 0)
        self.assertEqual([r["id_b"] for r in result.unmatched_b], ["b"])
        self.assertTrue(any("proper motion" in w
                            for w in result.warnings))

    def test_masked_coordinate_is_missing_not_phantom(self):
        from astropy.table import MaskedColumn
        a = Table()
        a["id"] = ["a"]
        # Underlying value sits exactly on the B row; the mask says the
        # measurement does not exist, so it must not match.
        a["ra"] = MaskedColumn([10.0], mask=[True])
        a["dec"] = [30.0]
        a["pmra"] = [0.0]
        a["pmdec"] = [0.0]
        result = match_tables(a, self._one(), epoch_a=2000.0,
                              epoch_b=2015.5, to_epoch=2016.0,
                              radius_arcsec=2.0, id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 0)
        self.assertEqual(len(result.unmatched_a), 1)
        self.assertTrue(any("non-finite" in w for w in result.warnings))

    def test_masked_pm_is_missing(self):
        from astropy.table import MaskedColumn
        b = Table()
        b["id"] = ["b"]
        b["ra"] = [10.0]
        b["dec"] = [30.0]
        b["pmra"] = MaskedColumn([50000.0], mask=[True])
        b["pmdec"] = [0.0]
        result = match_tables(self._one(), b, epoch_a=2000.0,
                              epoch_b=2015.5, to_epoch=2016.0,
                              radius_arcsec=2.0, id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 0)
        self.assertEqual(len(result.unmatched_b), 1)

    def test_nonfinite_radius_and_epochs_raise_clean_errors(self):
        from epochmatch.core import EpochMatchError
        for bad in (float("nan"), float("inf"), -1.0):
            with self.assertRaises(EpochMatchError):
                match_tables(self._one(), self._one(), epoch_a=2000.0,
                             epoch_b=2015.5, to_epoch=2016.0,
                             radius_arcsec=bad, id_a="id", id_b="id")
        with self.assertRaises(EpochMatchError):
            match_tables(self._one(), self._one(), epoch_a=float("nan"),
                         epoch_b=2015.5, to_epoch=2016.0,
                         radius_arcsec=2.0, id_a="id", id_b="id")

    def test_radius_above_180_degrees_still_matches(self):
        a = self._one()
        b = make_table(["y"], [160.0], [30.0], pmra=[0.0], pmdec=[0.0])
        # True separation is ~113.6 degrees; a 270-degree radius must
        # include it (the KD-tree chord must not wrap around).
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2000.0,
                              to_epoch=2000.0, radius_arcsec=270 * 3600.0,
                              id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 1)

    def test_pair_exactly_at_radius_matches(self):
        a = self._one()
        b = make_table(["y"], [10.0], [30.0 + 2.0 / 3600.0],
                       pmra=[0.0], pmdec=[0.0])
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2000.0,
                              to_epoch=2000.0, radius_arcsec=2.0,
                              id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 1)

    def test_b_row_cannot_be_matched_and_ambiguous_at_once(self):
        a = make_table(["a1", "a2"], [10.0, 10.0 + 1.5 / 3600.0],
                       [30.0, 30.0], pmra=[0.0, 0.0], pmdec=[0.0, 0.0])
        b = make_table(["b1", "b2"], [10.0, 10.0 - 1.0 / 3600.0],
                       [30.0, 30.0], pmra=[0.0, 0.0], pmdec=[0.0, 0.0])
        result = match_tables(a, b, epoch_a=2000.0, epoch_b=2000.0,
                              to_epoch=2000.0, radius_arcsec=2.0,
                              id_a="id", id_b="id")
        matched_b = {r["id_b"] for r in result.matched}
        amb_b = {r["id_b"] for r in result.ambiguous}
        self.assertFalse(matched_b & amb_b)
        self.assertEqual(len(result.matched), 0)
        self.assertIn("b1", amb_b)

    def test_dec_out_of_range_is_a_clean_error(self):
        from epochmatch.core import EpochMatchError
        a = make_table(["a"], [10.0], [95.0], pmra=[0.0], pmdec=[0.0])
        with self.assertRaises(EpochMatchError):
            match_tables(a, self._one(), epoch_a=2000.0, epoch_b=2015.5,
                         to_epoch=2016.0, radius_arcsec=2.0,
                         id_a="id", id_b="id")

    def test_partial_pm_columns_warn_which_was_found(self):
        b = make_table(["b"], [10.0], [30.0], pmra=[0.0])  # no pmdec
        result = match_tables(self._one(), b, epoch_a=2000.0,
                              epoch_b=2015.5, to_epoch=2016.0,
                              radius_arcsec=2.0, id_a="id", id_b="id")
        self.assertEqual(len(result.matched), 1)
        self.assertTrue(any("only 'pmra'" in w
                            for w in result.warnings))


if __name__ == "__main__":
    unittest.main()
