"""Command-line interface for epochmatch."""

from __future__ import annotations

import argparse
import sys

from .core import EpochMatchError, load_catalogue, match_tables


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="epochmatch",
        description="Cross-match two catalogues from different epochs. "
                    "Positions are propagated to a common epoch with "
                    "proper motions before matching, so moving stars "
                    "still match.")
    p.add_argument("catalogue_a", help="first catalogue (CSV, VOTable or "
                                       "FITS table)")
    p.add_argument("catalogue_b", help="second catalogue (CSV, VOTable or "
                                       "FITS table)")
    p.add_argument("--epoch-a", type=float, required=True,
                   help="epoch of catalogue A, as a year (e.g. 2000.0)")
    p.add_argument("--epoch-b", type=float, required=True,
                   help="epoch of catalogue B, as a year (e.g. 2015.5)")
    p.add_argument("--to-epoch", type=float, default=None,
                   help="epoch to propagate both catalogues to "
                        "(default: the later of the two epochs)")
    p.add_argument("--radius", type=float, default=2.0,
                   help="match radius in arcsec (default: 2.0)")
    p.add_argument("--out", default="epochmatch_result",
                   help="output prefix; writes <out>_matched.csv, "
                        "<out>_unmatched_a.csv, <out>_unmatched_b.csv and "
                        "<out>_ambiguous.csv (default: epochmatch_result)")
    p.add_argument("--ra-col", default="ra", help="RA column (degrees)")
    p.add_argument("--dec-col", default="dec", help="Dec column (degrees)")
    p.add_argument("--pmra-col", default="pmra",
                   help="proper-motion column in RA, mas/yr "
                        "(pm_ra_cosdec convention)")
    p.add_argument("--pmdec-col", default="pmdec",
                   help="proper-motion column in Dec, mas/yr")
    p.add_argument("--id-a", default=None,
                   help="ID column in catalogue A (default: row number)")
    p.add_argument("--id-b", default=None,
                   help="ID column in catalogue B (default: row number)")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    to_epoch = args.to_epoch
    if to_epoch is None:
        to_epoch = max(args.epoch_a, args.epoch_b)
    try:
        table_a = load_catalogue(args.catalogue_a)
        table_b = load_catalogue(args.catalogue_b)
        result = match_tables(
            table_a, table_b,
            epoch_a=args.epoch_a, epoch_b=args.epoch_b,
            to_epoch=to_epoch, radius_arcsec=args.radius,
            ra_col=args.ra_col, dec_col=args.dec_col,
            pmra_col=args.pmra_col, pmdec_col=args.pmdec_col,
            id_a=args.id_a, id_b=args.id_b)
        outputs = {
            "matched": result.matched,
            "unmatched_a": result.unmatched_a,
            "unmatched_b": result.unmatched_b,
            "ambiguous": result.ambiguous,
        }
        for kind, table in outputs.items():
            table.write(f"{args.out}_{kind}.csv", format="ascii.csv",
                        overwrite=True)
    except EpochMatchError as exc:
        print(f"epochmatch: error: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # never leak a traceback to the user
        print(f"epochmatch: error: {' '.join(str(exc).split())}",
              file=sys.stderr)
        return 2

    for warning in result.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    print(f"Matched: {len(result.matched)}")
    print(f"Ambiguous: {len(result.ambiguous)}")
    print(f"Unmatched in A: {len(result.unmatched_a)}")
    print(f"Unmatched in B: {len(result.unmatched_b)}")
    print(f"Wrote: {args.out}_matched.csv, {args.out}_unmatched_a.csv, "
          f"{args.out}_unmatched_b.csv, {args.out}_ambiguous.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
