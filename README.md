# epochmatch — cross-match star catalogues from different epochs

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/RakTikoo/epochmatch)

**Cross-match two star catalogues from different epochs — positions are moved to a common epoch with proper motions first, so stars that moved still match.**

Matching two catalogues by raw position only works if both were measured at the same time. They usually were not: a star with a proper motion of 400 mas/yr sits ~6 arcsec away from its year-2000 position by 2016, far outside a normal 1–2 arcsec match radius. The usual matchers treat the sky as frozen, so movers silently mismatch or vanish. epochmatch propagates both catalogues to one shared epoch before matching, and tells you plainly which pairs matched, which did not, and which were too close to call.

**Quick install:** `pip install git+https://github.com/RakTikoo/epochmatch` (Python 3.10+; Astropy, NumPy and SciPy). Details under [Install](#install).

## Demo

`examples/cat_a.csv` holds 6 stars at epoch 2000.0. `examples/cat_b.fits` holds the same field at epoch 2015.5: five of the stars (moved by their proper motions), one star missing, one extra star, and two extra stars planted close to S3.

```console
$ epochmatch examples/cat_a.csv examples/cat_b.fits \
    --epoch-a 2000.0 --epoch-b 2015.5 --to-epoch 2016.0 \
    --radius 2.0 --id-a id --id-b id --out demo_result
Matched: 4
Ambiguous: 3
Unmatched in A: 1
Unmatched in B: 1
Wrote: demo_result_matched.csv, demo_result_unmatched_a.csv, demo_result_unmatched_b.csv, demo_result_ambiguous.csv
```

`demo_result_matched.csv` (S3 is not here — see below):

```csv
id_a,id_b,sep_arcsec,a_ra,a_dec,a_pmra,a_pmdec,b_ra,b_dec,b_pmra,b_pmdec
S1,S1,0.0,150.0,20.0,120.0,-50.0,150.00054982443052,19.99978472137435,120.0,-50.0
S2,S2,0.0,150.05,20.04,-80.0,90.0,150.04963335572077,20.040387499622337,-80.0,90.0
S4,S4,0.0,150.15,20.12,40.0,10.0,150.15018341533047,20.120043055460734,40.0,10.0
S5,S5,0.0,150.2,20.16,-200.0,150.0,150.19908268497232,20.16064583095752,-200.0,150.0
```

S3 moved ~6 arcsec between the epochs, and at the common epoch it has three B rows within the radius (its true counterpart plus the two planted neighbours), so all three pairs land in `demo_result_ambiguous.csv` with the reason `multiple_candidates_within_radius` instead of a guessed match. S6 is in `demo_result_unmatched_a.csv`; the extra star X1 is in `demo_result_unmatched_b.csv`.

## Install

```console
$ pip install git+https://github.com/RakTikoo/epochmatch
```

Needs Python 3.10+; depends on Astropy, NumPy and SciPy. Or from a clone: `pip install .`

Or install the wheel straight from this repo's release assets (no package
registry involved):

```sh
pip install https://github.com/RakTikoo/epochmatch/releases/download/v0.1.0/epochmatch-0.1.0-py3-none-any.whl
```

Or install from the git repo with pipx or uv:

```sh
pipx install git+https://github.com/RakTikoo/epochmatch.git@v0.1.0
uv tool install --from git+https://github.com/RakTikoo/epochmatch.git@v0.1.0 epochmatch
```

With Homebrew (macOS/Linux):

```sh
brew tap RakTikoo/tools
brew install RakTikoo/tools/epochmatch
```

With Scoop (Windows):

```powershell
scoop bucket add RakTikoo https://github.com/RakTikoo/scoop-bucket
scoop install epochmatch
```

**For agent users:** this repository ships a `SKILL.md` for AI agents. Install it with `npx skills add RakTikoo/epochmatch`.

## Usage

```console
$ epochmatch A.csv B.fits --epoch-a 2000.0 --epoch-b 2015.5 --to-epoch 2016.0 --radius 2.0 --out result
```

- Inputs can be CSV, VOTable or FITS tables (anything `astropy.table` reads), mixed freely.
- Columns: `ra`, `dec` in degrees; `pmra`, `pmdec` in mas/yr, using the `pm_ra_cosdec` convention (as in Gaia). Rename any of them with `--ra-col`, `--dec-col`, `--pmra-col`, `--pmdec-col`.
- `--id-a` / `--id-b` name an ID column per catalogue. Without one, the ID is the row number (0-based).
- `--propagated-columns` adds the propagated position of each side at the common epoch to the matched output (`ra_a_prop_deg`, `dec_a_prop_deg`, `ra_b_prop_deg`, `dec_b_prop_deg`). Default behaviour is unchanged when the flag is omitted.
- `--to-epoch` defaults to the later of the two epochs.
- If a catalogue has no proper-motion columns, its positions are matched as-is and a clear warning is printed — the match for that side is static, and the summary never hides that.
- Rows with missing or non-finite coordinates or proper motions cannot match; they are counted in a warning and reported as unmatched. Masked/null values (as in Gaia-style FITS and VOTable files) count as missing — the mask is honoured, never silently dropped.

### How a pair is decided

Both catalogues are propagated to the common epoch (Astropy `apply_space_motion`). For each A row, all B rows inside `--radius` are candidates:

- exactly one candidate → **matched**, unless another A row claims the same B row, in which case all of those pairs are **ambiguous**;
- two or more candidates → **ambiguous** (every candidate pair is listed, with the reason);
- no candidate → **unmatched**.

### Outputs

All CSV, written next to wherever you run the command:

| File | Contents |
|---|---|
| `<out>_matched.csv` | `id_a`, `id_b`, `sep_arcsec` (at the common epoch), plus both rows' coordinate and proper-motion columns |
| `<out>_unmatched_a.csv` / `<out>_unmatched_b.csv` | rows with no candidate, with their coordinates |
| `<out>_ambiguous.csv` | one row per candidate pair, with separation and reason |

## Limits

Honest list:

- **Directions only.** Propagation uses proper motions; parallax and radial velocity are not used (an assumed distance cancels out of the propagated direction when radial velocity is zero). For very nearby, very fast stars over very long baselines, a full space-motion treatment with measured parallaxes would be more accurate.
- **Units are fixed by convention:** degrees and mas/yr. There are no unit-conversion flags; convert your columns first if they differ.
- **The proper-motion convention is `pm_ra_cosdec`** (the RA motion already includes the cos(Dec) factor), matching Gaia and most modern catalogues.
- **Ambiguity is positional only.** epochmatch does not do probabilistic or magnitude-aware matching (no NWAY-style Bayesian cross-match); close pairs are reported, not resolved.
- **Not built for billion-row surveys.** Matching is a KD-tree over both catalogues in memory — comfortable for thousands to millions of rows, not a Gaia-scale join. For that, use a partitioned tool such as LSDB.
- **Dense fields grow quadratically in the output.** If N stars on each side all fall inside one match radius, the ambiguous file gets N×N candidate rows — that is the honest answer, but it is big. (Measured: a 300×300 fully coincident field finishes in about a second; a sparse 20,000×20,000 match takes about 2 seconds.)
- **Extreme proper motions are model-sensitive.** Propagation is exact linear space motion (via Astropy/ERFA). A constant-angular-rate approximation differs from it by ~0.3 arcsec at |pm| ≈ 100,000 mas/yr over a 16-year baseline — but no known star moves that fast (Barnard's star, the fastest, is ~10,300 mas/yr, where the difference is ~0.003 arcsec), so this does not affect real catalogues.
- Separations are computed at the common epoch; by default the coordinate columns in the outputs are the original catalogue values, not the propagated ones. Pass `--propagated-columns` to also write the propagated positions (at the common epoch) as `ra_a_prop_deg`, `dec_a_prop_deg`, `ra_b_prop_deg` and `dec_b_prop_deg` in the matched output. For a catalogue without proper-motion columns, the proper-motion columns in the matched output are NaN rather than a fabricated 0.0.

## Tests

```console
$ python3 run_tests.py
```

22 tests, including: proper-motion stars a static match provably misses but epochmatch recovers; zero-motion agreement with a plain `SkyCoord` match; ambiguous pairs; missing proper-motion columns; empty catalogues; NaN coordinates; duplicate positions; radius 0; CLI error handling (one clean line, no traceback); and the adversarial regressions from `REVIEW.md` — NaN/masked proper motions, masked coordinates, non-finite radius and epochs, a pair exactly at the radius, radii above 180 degrees, a B row claimed by both an ambiguous and a single-candidate A row, and out-of-range Dec.

## Use with an AI agent

This repo ships a [SKILL.md](SKILL.md) written to the Agent Skills open
standard. If your coding agent supports skills, point it at this repo and
it can install epochmatch, run it on your catalogues, and read the results
for you. The skill works as documented when the repo is cloned under its
real name (`epochmatch`) — the standard matches a skill to its folder name.

## License

MIT — see `LICENSE`.
