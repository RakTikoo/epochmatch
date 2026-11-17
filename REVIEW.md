# REVIEW.md — independent adversarial review of epochmatch

Reviewer: independent subagent (did not build the tool). Date: 2026-11-17.
Scope: local work only; nothing published or pushed. All fixes are inside
this directory (`epochmatch/core.py`, `tests/test_epochmatch.py`,
`README.md`, this file).

## Headline-claim verification (independent oracle)

The claim: propagating both catalogues to a common epoch recovers matches
that a naive position-only match misses.

Independent check, written from first principles with numpy trig only
(no epochmatch code, no `apply_space_motion`):

- Synthetic star (RA 10°, Dec 30°, pmra 500, pmdec 300 mas/yr, epoch
  2000.0): manual great-circle propagation puts it 9.038 arcsec away by
  2015.5. Naive separation A(2000) vs B(2015.5) = 9.038 arcsec → a
  2-arcsec naive match FAILS. Both sides manually propagated to 2016.0:
  separation 0.000006 arcsec → epoch match SUCCEEDS. epochmatch on the
  same fixture: matched, sep 0.0.
- Shipped demo files: the oracle reproduces every `cat_b.fits` position
  from `cat_a.csv` with 0.000 mas residual (S1, S3, S5 checked). S3's
  displacement is 6.053 arcsec, confirming the README's "~6 arcsec".
- High-PM check with an exact *linear* space-motion oracle (constant
  velocity vector; B's position and evolved proper motion at 2015.5
  derived independently): epochmatch matches at pm = 400, 10,300
  (Barnard's star), 100,000 and 1,000,000 mas/yr with separations
  0.0, 0.0, 0.0002 and 0.1424 arcsec. See F12 for why the oracle model
  matters.

Verdict: headline claim TRUE and independently verified.

## Fresh install + demo

- Fresh venv (`python3 -m venv`, `pip install ~/workspace/astro-tools/epochmatch`):
  installs cleanly (astropy 8.0.2, numpy 2.5.5, scipy 1.18.1), console
  script works. Done twice — before and after fixes.
- README demo reproduced EXACTLY, before and after fixes: stdout
  `Matched: 4 / Ambiguous: 3 / Unmatched in A: 1 / Unmatched in B: 1`,
  and all four output CSVs byte-identical to the pre-fix run, with the
  matched CSV identical to the CSV printed in the README.

## Tests

- Original suite: 12 tests, all passing (`run_tests.py`, astro-tools venv).
- After fixes: 22 tests, all passing (12 original + 10 new regression
  tests in `TestAdversarialRegressions`), verified with both the
  astro-tools venv python and a fresh-venv python.

## Findings and fixes

### F1 — HIGH — Masked/null values were silently unmasked
`np.asarray(column, dtype=float)` drops astropy masks and exposes the
underlying stored value. A masked RA whose underlying value sat exactly
on a B row produced a phantom match with no warning; a masked PM
propagated its underlying garbage value (50000 mas/yr in the test). A
real FITS file with a TNULL-masked PM (Gaia-style) crashed the whole run
(ERFA non-convergence, then a scipy "must be finite" error) — one null
row destroyed every match.
Fix: `_float_column()` converts masked entries to NaN, so they flow into
the missing-data path (warning + unmatched). Regression tests added.

### F2 — HIGH — One row with NaN/Inf proper motion crashed everything
Finiteness was checked for RA/Dec only. A single NaN PM in A raised
`ValueError: 'x' must be finite` from scipy; in B, `ValueError: data
must be finite` from the KD-tree. Via the API these were raw tracebacks;
via the CLI the whole run was lost to one cryptic line. An Inf PM was
worse: ERFA clamped it and the row silently matched a static position.
Fix: rows with non-finite PMs (when PM columns exist) are treated
exactly like bad coordinates — excluded, counted in a warning
("…have missing or non-finite proper motions; they cannot be
propagated and are reported as unmatched"). Regression tests added.

### F3 — MEDIUM — NaN/Inf radius silently "succeeded" with zero matches
`--radius nan` and `--radius inf` exited 0 printing `Matched: 0` —
a silently wrong success. (Negative radius was already a clean error.)
Fix: `match_tables` validates the radius is finite and >= 0, raising
`EpochMatchError`. CLI now exits 2 with
`epochmatch: error: radius must be finite and >= 0 arcsec`.

### F4 — MEDIUM — A B row could be matched and ambiguous at the same time
With A1 seeing {B1, B2} (ambiguous) and A2 seeing only {B1}, the old
code matched A2–B1 while also listing B1 as ambiguous — contradicting
the README's own rule ("exactly one candidate → matched, unless another
A row claims the same B row") and producing self-contradictory outputs.
Fix: a tentative pair whose B row is already a candidate of an
ambiguous A row becomes ambiguous with reason
`b_row_claimed_by_multiple_a_rows`. Demo output unchanged. Test added.

### F5 — MEDIUM — Performance made the README scale claim untrue
Separations were computed with one scalar `SkyCoord.separation` call
per candidate pair (~2 ms each). Measured before the fix: sparse
20,000×20,000 identical catalogues took 46.0 s; a dense 300×300 fully
coincident field took 181.8 s. "Comfortable for thousands to millions
of rows" was not credible.
Fix: separations are computed vectorised from the chord length,
sep = 2·asin(chord/2) (well conditioned at small angles, unlike
arccos of the dot product). After: 1.64 s and 0.98 s respectively,
with demo outputs byte-identical. README Limits now also states the
inherent dense-field caveat (N×N ambiguous rows are the output itself).

### F6 — LOW — Radii above 180° missed real matches
The KD-tree chord `2·sin(θ/2)` wraps for θ > 360° of radius: with a
270° radius, a pair 113.6° apart was not matched.
Fix: the query radius is capped at the whole sphere (the exact
separation filter still enforces the requested radius). Test added.

### F7 — LOW — A pair exactly at the radius was dropped
Floating-point rounding in the cartesian prefilter dropped a pair at
exactly 2.000000 arcsec with radius 2.0.
Fix: the prefilter queries with a padded radius (+1 µas + 1e-12
relative); the exact filter is unchanged. Test added.

### F8 — LOW — Non-finite epochs leaked raw astropy errors via the API
`epoch_a=nan` / `to_epoch=inf` raised a multi-line `ValueError` from
`astropy.time`, not an `EpochMatchError` (the CLI caught it, but with a
cryptic message). Fix: epochs are validated in `match_tables` (and in
`propagate_positions`); CLI now prints e.g.
`epochmatch: error: epoch_a must be a finite year, got nan`. Test added.

### F9 — LOW — Dec outside [-90, +90] leaked a raw astropy error
Dec = 95° raised `ValueError: Latitude angle(s) must be within…` from
SkyCoord construction. Fix: `Catalogue` raises `EpochMatchError`
("Dec column … has values outside -90 to +90 degrees"). Test added.
(RA is deliberately not range-checked: any finite RA is a well-defined
direction modulo 360°, and matching stays correct — verified.)

### F10 — LOW — Fabricated 0.0 proper motions in the matched output
For a catalogue with no PM columns, the matched CSV reported
`b_pmra = 0.0` — a measurement the catalogue never made.
Fix: those columns are now NaN for a no-PM side; README says so.

### F11 — LOW — Misleading warning when only one PM column exists
With only `pmra` present, the warning claimed the catalogue "has no
proper-motion columns ('pmra'/'pmdec')". Fix: the warning now appends
"(only 'pmra' was found; both columns are required, so it was not
used)". Behaviour (static match + warning) was already safe. Test added.

### F12 — INVESTIGATED, NOT A BUG — Extreme proper motions
Suspicion: at pm = 10⁶ mas/yr epochmatch disagreed with a great-circle
oracle by 58 arcsec. Investigation: the *oracle* was the approximation.
Astropy/ERFA `apply_space_motion` performs exact linear space motion
(constant velocity vector); a constant-angular-rate great circle differs
from it by ≈ x³/3 (x = total motion in radians) — 0.003 arcsec at
Barnard's-star rates over 16 yr, 0.3 arcsec at 10⁵ mas/yr, 58 arcsec at
10⁶ mas/yr. With an exact linear oracle, epochmatch is correct at all
tested rates (see headline section). No real star is affected (fastest
known ≈ 10,300 mas/yr). Documented in README Limits.
Related: ERFA's "distance overridden" note fires for high-PM stars at
the assumed 1 kpc distance. It is benign here — with zero radial
velocity the assumed distance cancels out of the direction (verified:
identical propagated directions at 1 pc, 10 pc, 100 pc, 1 kpc) — so
`Catalogue.propagated` now filters that one note out, while genuine
ERFA warnings (non-convergence, excessive velocity) are collected and
surfaced as tool warnings instead of only appearing as Python warnings.

## Checked and found correct (no change)

- RA 0h seam pair (359.9999° vs 0.0001°) matches at 0.72 arcsec;
  near-pole positions (Dec 89.9999°) with PM match.
- Zero-motion agreement with a plain SkyCoord match (existing test).
- Both catalogues empty; header-only CSV; VOTable load; string PM and
  sexagesimal RA columns → clean one-line errors; missing file, missing
  column, directory-as-input, unwritable output path → exit 2, one line,
  no traceback (CLI).
- Duplicate positions → ambiguous; radius 0 → exact coincidences only.
- README's quantitative claims: "400 mas/yr … ~6 arcsec by 2016"
  (computed 6.4"), S3 "~6 arcsec" (6.053"), "an assumed distance cancels
  out of the propagated direction when radial velocity is zero"
  (verified numerically across 1 pc–1 kpc).

## README assessment

The Demo, Usage and decision-rule sections were accurate and reproduced
exactly. The Limits section was honest but incomplete in three places,
now amended: (1) the missing-data sentence covered coordinates only —
it now covers proper motions and masked/null values; (2) no mention of
dense-field quadratic ambiguous output — added, with measured timings;
(3) no mention of extreme-PM model sensitivity — added, quantified.
Test count updated 12 → 22. One claim could not be verified under the
local-only constraint of this review: the install line
`pip install git+https://github.com/RakTikoo/epochmatch` (the local
equivalent, `pip install .` from a fresh venv, was verified twice).

## Not fixed (by design or out of scope)

- Duplicate ID values in a catalogue make output rows indistinguishable
  — that is the caller's data, not a matching error.
- Sexagesimal/string coordinates and non-degree units are rejected, as
  documented ("Units are fixed by convention").
- No probabilistic/magnitude-aware matching — documented limit, a
  design decision, left unchanged.
