#!/usr/bin/env python
import sys
import os
import argparse
import re
from astropy.coordinates import SkyCoord
from astropy import units as u
import pulsarsurveyscraper


def main():
    pulsarsurveyscraper.log.setLevel("WARNING")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("coord", nargs="+", help="Coordinates to search")
    parser.add_argument("-d", "--dm", type=float, help="DM to match")
    parser.add_argument("--dmtol", default=10, type=float, help="DM tolerance")
    parser.add_argument(
        "-r", "--radius", default=5, type=float, help="Search radius (degrees)"
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        dest="dest",
        default="",
        help="Input directory for survey caches",
    )
    parser.add_argument(
        "-v", "--verbosity", default=0, action="count", help="Increase output verbosity"
    )

    args = parser.parse_args()
    if args.verbosity == 1:
        pulsarsurveyscraper.log.setLevel("INFO")
    elif args.verbosity >= 2:
        pulsarsurveyscraper.log.setLevel("DEBUG")

    if len(args.coord) == 2:
        ra, dec = args.coord
    elif len(args.coord) == 1:
        c = args.coord[0].split()
        l = len(c)
        ra = " ".join(c[: (l // 2)])
        dec = " ".join(c[(l // 2) :])

    try:
        if (re.search(r"[^\d.+\-]", ra) is None) and (
            re.search(r"[^\d.+\-]", dec) is None
        ):
            coord = SkyCoord(ra, dec, unit="deg")
        else:
            coord = SkyCoord(ra, dec)
    except ValueError:
        try:
            coord = SkyCoord(ra, dec, unit=("hour", "deg"))
        except ValueError:
            pulsarsurveyscraper.log.error(
                "Unable to parse input coordinates '{},{}'".format(ra, dec)
            )
            sys.exit(1)

    pulsar_table = pulsarsurveyscraper.PulsarTable(directory=args.dest)
    print(
        "Searching {:.1f}deg around {} = {}d,{}d".format(
            args.radius,
            coord.to_string("hmsdms", sep=":"),
            coord.ra.to_string(decimal=True),
            coord.dec.to_string(decimal=True, alwayssign=True),
        )
    )
    if args.dm is not None and args.dmtol is not None:
        print(
            "Also requiring DM with +/-{:.1f} of {:.1f} pc/cm**2".format(
                args.dmtol, args.dm
            )
        )
    result = pulsar_table.search(
        coord, radius=args.radius * u.deg, DM=args.dm, DM_tolerance=args.dmtol
    )
    print("Found {} matches:".format(len(result)))
    print(result)


if __name__ == "__main__":
    main()
