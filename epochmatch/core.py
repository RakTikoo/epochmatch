"""Core matching logic for epochmatch.

Cross-match two catalogues whose positions were measured at different
epochs.  Both catalogues are propagated to a common epoch with their
proper-motion columns before matching, because stars move: matching raw
positions from different epochs mismatches anything that moved more than
the match radius between the two epochs.
"""

from __future__ import annotations

import os
import warnings as _warnings
from collections import Counter

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import MaskedColumn, Table
from astropy.time import Time
from scipy.spatial import cKDTree


class EpochMatchError(Exception):
    """A user-facing error.  Printed as one line, never a traceback."""


def _one_line(msg) -> str:
    return " ".join(str(msg).split())


def _float_column(table: Table, name: str) -> np.ndarray:
    """A table column as a plain float array.

    Masked entries (nulls in FITS/VOTable files, e.g. Gaia-style missing
    values) become NaN, so the caller's finiteness checks see them.
    ``np.asarray`` alone would silently drop the mask and expose the
    underlying stored value, inventing data the catalogue says is
    missing.
    """
    col = table[name]
    arr = np.asarray(col, dtype=float)
    if isinstance(col, MaskedColumn):
        mask = np.ma.getmaskarray(col)
        if mask.any():
            arr = arr.copy()
            arr[mask] = np.nan
    return arr


def load_catalogue(path: str) -> Table:
    if not os.path.isfile(path):
        raise EpochMatchError(f"file not found: {path}")
    try:
        return Table.read(path)
    except Exception as exc:  # astropy raises many different types here
        raise EpochMatchError(
            f"could not read {path} as a table (CSV, VOTable or FITS "
            f"expected): {_one_line(exc)}"
        )


class Catalogue:
    """One input catalogue, normalised to plain arrays."""

    def __init__(self, table: Table, label: str, ra_col: str, dec_col: str,
                 pmra_col: str, pmdec_col: str, id_col: str | None):
        self.label = label
        self.warnings: list[str] = []
        n = len(table)

        for name, role in ((ra_col, "RA"), (dec_col, "Dec")):
            if name not in table.colnames:
                raise EpochMatchError(
                    f"catalogue {label}: missing {role} column '{name}' "
                    f"(has: {', '.join(table.colnames)})"
                )
        if id_col is not None and id_col not in table.colnames:
            raise EpochMatchError(
                f"catalogue {label}: missing ID column '{id_col}' "
                f"(has: {', '.join(table.colnames)})"
            )

        try:
            self.ra = _float_column(table, ra_col)
            self.dec = _float_column(table, dec_col)
        except (TypeError, ValueError):
            raise EpochMatchError(
                f"catalogue {label}: RA/Dec columns must be numeric degrees"
            )

        bad_dec = np.isfinite(self.dec) & (np.abs(self.dec) > 90.0)
        if bad_dec.any():
            raise EpochMatchError(
                f"catalogue {label}: Dec column '{dec_col}' has values "
                "outside -90 to +90 degrees"
            )

        if pmra_col in table.colnames and pmdec_col in table.colnames:
            try:
                self.pmra = _float_column(table, pmra_col)
                self.pmdec = _float_column(table, pmdec_col)
            except (TypeError, ValueError):
                raise EpochMatchError(
                    f"catalogue {label}: proper-motion columns must be "
                    "numeric mas/yr"
                )
            self.has_pm = True
        else:
            self.pmra = np.zeros(n)
            self.pmdec = np.zeros(n)
            self.has_pm = False
            msg = (
                f"catalogue {label} has no proper-motion columns "
                f"('{pmra_col}'/'{pmdec_col}'); its positions were matched "
                "statically, with no epoch propagation"
            )
            present = [c for c in (pmra_col, pmdec_col)
                       if c in table.colnames]
            if len(present) == 1:
                msg += (f" (only '{present[0]}' was found; both columns "
                        "are required, so it was not used)")
            self.warnings.append(msg)

        if id_col is not None:
            self.ids = [str(v) for v in table[id_col]]
        else:
            self.ids = [str(i) for i in range(n)]

        # Rows with non-finite coordinates can never be matched.  Keep them
        # out of the matching, count them, and let them land in unmatched.
        self.finite = np.isfinite(self.ra) & np.isfinite(self.dec)
        n_bad = int((~self.finite).sum())
        if n_bad:
            self.warnings.append(
                f"catalogue {label}: {n_bad} row(s) have missing or "
                "non-finite coordinates; they cannot match and are "
                "reported as unmatched"
            )

        # Rows with proper-motion columns present but non-finite values
        # (NaN, inf, or masked/null) cannot be propagated: propagating
        # them either crashes the matcher or, worse, silently invents a
        # position.  Treat them exactly like bad coordinates.
        if self.has_pm:
            pm_finite = np.isfinite(self.pmra) & np.isfinite(self.pmdec)
            n_pm_bad = int((~pm_finite & self.finite).sum())
            if n_pm_bad:
                self.warnings.append(
                    f"catalogue {label}: {n_pm_bad} row(s) have missing "
                    "or non-finite proper motions; they cannot be "
                    "propagated and are reported as unmatched"
                )
            self.finite = self.finite & pm_finite

    def propagated(self, from_epoch: float, to_epoch: float) -> SkyCoord:
        """Positions of the finite rows, propagated to ``to_epoch``."""
        idx = np.nonzero(self.finite)[0]
        coord = SkyCoord(
            ra=self.ra[idx] * u.deg,
            dec=self.dec[idx] * u.deg,
            pm_ra_cosdec=self.pmra[idx] * u.mas / u.yr,
            pm_dec=self.pmdec[idx] * u.mas / u.yr,
            distance=1.0 * u.kpc,
            radial_velocity=0.0 * u.km / u.s,
            obstime=Time(from_epoch, format="jyear"),
        )
        if to_epoch != from_epoch and len(idx) > 0:
            with _warnings.catch_warnings(record=True) as caught:
                _warnings.simplefilter("always")
                coord = coord.apply_space_motion(
                    new_obstime=Time(to_epoch, format="jyear"))
            seen: list[str] = []
            for w in caught:
                msg = _one_line(w.message)
                # ERFA's "distance overridden" note is benign here by
                # construction: with zero radial velocity the assumed
                # distance cancels out of the propagated direction
                # (verified against an independent calculation), so
                # ERFA substituting its own distance changes nothing.
                # Other warnings (non-convergence, excessive velocity)
                # are real and are surfaced below.
                if "distance overridden" in msg:
                    continue
                if msg not in seen:
                    seen.append(msg)
            for msg in seen[:3]:
                self.warnings.append(
                    f"catalogue {self.label}: astropy/ERFA propagation "
                    f"warning: {msg}"
                )
        return coord


def propagate_positions(ra, dec, pmra, pmdec, from_epoch, to_epoch):
    """Propagate plain arrays to a new epoch; returns a SkyCoord.

    Used by the matcher and handy in tests.  The assumed distance is
    arbitrary: with zero radial velocity the propagated *direction*
    does not depend on it.
    """
    for name, val in (("from_epoch", from_epoch), ("to_epoch", to_epoch)):
        if not np.isfinite(val):
            raise EpochMatchError(f"{name} must be a finite year, "
                                  f"got {val}")
    coord = SkyCoord(
        ra=np.asarray(ra, dtype=float) * u.deg,
        dec=np.asarray(dec, dtype=float) * u.deg,
        pm_ra_cosdec=np.asarray(pmra, dtype=float) * u.mas / u.yr,
        pm_dec=np.asarray(pmdec, dtype=float) * u.mas / u.yr,
        distance=1.0 * u.kpc,
        radial_velocity=0.0 * u.km / u.s,
        obstime=Time(from_epoch, format="jyear"),
    )
    if to_epoch != from_epoch:
        coord = coord.apply_space_motion(
            new_obstime=Time(to_epoch, format="jyear"))
    return coord


def _candidates(prop_a: SkyCoord, prop_b: SkyCoord, radius_arcsec: float):
    """For each A row, the B rows within radius, as (b_index, sep_arcsec)."""
    empty = [[] for _ in range(len(prop_a))]
    if len(prop_a) == 0 or len(prop_b) == 0:
        return empty
    va = np.vstack([prop_a.cartesian.x.value,
                    prop_a.cartesian.y.value,
                    prop_a.cartesian.z.value]).T
    vb = np.vstack([prop_b.cartesian.x.value,
                    prop_b.cartesian.y.value,
                    prop_b.cartesian.z.value]).T
    va /= np.linalg.norm(va, axis=1, keepdims=True)
    vb /= np.linalg.norm(vb, axis=1, keepdims=True)
    tree = cKDTree(vb)
    # Query with a slightly padded radius, capped at the whole sphere:
    # floating-point rounding in the cartesian vectors can otherwise
    # drop a pair sitting exactly on the radius, and an uncapped chord
    # wraps around for radii above 180 degrees.  The exact separation
    # filter below still enforces the requested radius.
    query_radius = radius_arcsec + 1e-6 + abs(radius_arcsec) * 1e-12
    theta = min(np.radians(query_radius / 3600.0), np.pi)
    chord = 2.0 * np.sin(theta / 2.0)
    result = []
    for i, cand in enumerate(tree.query_ball_point(va, r=chord)):
        lst = []
        if cand:
            cand_arr = np.asarray(cand)
            # Separation from the chord length, sep = 2*asin(chord/2):
            # vectorised, and well conditioned at small separations
            # (unlike arccos of the dot product).
            diffs = vb[cand_arr] - va[i]
            chords = np.linalg.norm(diffs, axis=1)
            seps = np.degrees(
                2.0 * np.arcsin(np.clip(chords / 2.0, 0.0, 1.0))) * 3600.0
            for j, sep in zip(cand_arr, seps):
                if sep <= radius_arcsec + 1e-9:
                    lst.append((int(j), float(sep)))
        lst.sort(key=lambda t: t[1])
        result.append(lst)
    return result


class MatchResult:
    def __init__(self, matched, unmatched_a, unmatched_b, ambiguous,
                 warnings):
        self.matched = matched
        self.unmatched_a = unmatched_a
        self.unmatched_b = unmatched_b
        self.ambiguous = ambiguous
        self.warnings = warnings


def _key_columns(cat: Catalogue, rows, prefix: str) -> dict:
    if cat.has_pm:
        pmra, pmdec = cat.pmra[rows], cat.pmdec[rows]
    else:
        # No proper-motion columns in the source: report NaN rather
        # than a fabricated 0.0 the catalogue never measured.
        pmra = np.full(len(rows), np.nan)
        pmdec = np.full(len(rows), np.nan)
    return {
        f"{prefix}_ra": cat.ra[rows],
        f"{prefix}_dec": cat.dec[rows],
        f"{prefix}_pmra": pmra,
        f"{prefix}_pmdec": pmdec,
    }


def match_tables(table_a: Table, table_b: Table, *,
                 epoch_a: float, epoch_b: float, to_epoch: float,
                 radius_arcsec: float,
                 ra_col="ra", dec_col="dec", pmra_col="pmra",
                 pmdec_col="pmdec", id_a=None, id_b=None) -> MatchResult:
    try:
        radius_arcsec = float(radius_arcsec)
    except (TypeError, ValueError):
        raise EpochMatchError("radius must be a number of arcsec")
    if not np.isfinite(radius_arcsec) or radius_arcsec < 0:
        raise EpochMatchError("radius must be finite and >= 0 arcsec")
    for name, val in (("epoch_a", epoch_a), ("epoch_b", epoch_b),
                      ("to_epoch", to_epoch)):
        try:
            finite_epoch = np.isfinite(float(val))
        except (TypeError, ValueError):
            finite_epoch = False
        if not finite_epoch:
            raise EpochMatchError(f"{name} must be a finite year, "
                                  f"got {val}")
    cat_a = Catalogue(table_a, "A", ra_col, dec_col, pmra_col, pmdec_col,
                      id_a)
    cat_b = Catalogue(table_b, "B", ra_col, dec_col, pmra_col, pmdec_col,
                      id_b)

    prop_a = cat_a.propagated(epoch_a, to_epoch)
    prop_b = cat_b.propagated(epoch_b, to_epoch)
    # Collected after propagation: propagating can add warnings.
    warnings = cat_a.warnings + cat_b.warnings
    a_rows = np.nonzero(cat_a.finite)[0]   # finite-row index -> table row
    b_rows = np.nonzero(cat_b.finite)[0]

    cand = _candidates(prop_a, prop_b, radius_arcsec)

    tentative = {}
    ambiguous = []  # (a_row, b_row, sep, reason)
    amb_b_prop = set()  # B rows (propagated indexing) already ambiguous
    for i, lst in enumerate(cand):
        if len(lst) >= 2:
            for j, sep in lst:
                ambiguous.append(
                    (int(a_rows[i]), int(b_rows[j]), sep,
                     "multiple_candidates_within_radius"))
                amb_b_prop.add(j)
        elif len(lst) == 1:
            tentative[i] = lst[0]

    claims = Counter(j for j, _ in tentative.values())
    matched = []  # (a_row, b_row, sep)
    for i, (j, sep) in tentative.items():
        # A B row that is also a candidate of an ambiguous A row is
        # claimed by more than one A row in every sense that matters:
        # it cannot be a clean match here while listed as ambiguous
        # there (the README rule: exactly one candidate -> matched,
        # unless another A row claims the same B row).
        if claims[j] > 1 or j in amb_b_prop:
            ambiguous.append(
                (int(a_rows[i]), int(b_rows[j]), sep,
                 "b_row_claimed_by_multiple_a_rows"))
        else:
            matched.append((int(a_rows[i]), int(b_rows[j]), sep))

    matched_a = {a for a, _, _ in matched}
    matched_b = {b for _, b, _ in matched}
    amb_a = {a for a, _, _, _ in ambiguous}
    amb_b = {b for _, b, _, _ in ambiguous}
    unmatched_a_rows = [r for r in range(len(cat_a.ids))
                        if r not in matched_a and r not in amb_a]
    unmatched_b_rows = [r for r in range(len(cat_b.ids))
                        if r not in matched_b and r not in amb_b]

    matched_t = Table()
    if matched:
        ma = [a for a, _, _ in matched]
        mb = [b for _, b, _ in matched]
        matched_t["id_a"] = [cat_a.ids[a] for a in ma]
        matched_t["id_b"] = [cat_b.ids[b] for b in mb]
        matched_t["sep_arcsec"] = [round(s, 4) for _, _, s in matched]
        for k, v in _key_columns(cat_a, ma, "a").items():
            matched_t[k] = v
        for k, v in _key_columns(cat_b, mb, "b").items():
            matched_t[k] = v
    else:
        matched_t = Table(names=["id_a", "id_b", "sep_arcsec",
                                 "a_ra", "a_dec", "a_pmra", "a_pmdec",
                                 "b_ra", "b_dec", "b_pmra", "b_pmdec"])

    ambiguous_t = Table()
    if ambiguous:
        ambiguous_t["id_a"] = [cat_a.ids[a] for a, _, _, _ in ambiguous]
        ambiguous_t["id_b"] = [cat_b.ids[b] for _, b, _, _ in ambiguous]
        ambiguous_t["sep_arcsec"] = [round(s, 4) for _, _, s, _ in ambiguous]
        ambiguous_t["reason"] = [r for _, _, _, r in ambiguous]
    else:
        ambiguous_t = Table(names=["id_a", "id_b", "sep_arcsec", "reason"])

    def _unmatched(cat, rows, id_name, prefix):
        t = Table()
        t[id_name] = [cat.ids[r] for r in rows]
        t[f"{prefix}_ra"] = cat.ra[rows] if rows else np.array([])
        t[f"{prefix}_dec"] = cat.dec[rows] if rows else np.array([])
        return t

    return MatchResult(
        matched_t,
        _unmatched(cat_a, unmatched_a_rows, "id_a", "a"),
        _unmatched(cat_b, unmatched_b_rows, "id_b", "b"),
        ambiguous_t,
        warnings,
    )
