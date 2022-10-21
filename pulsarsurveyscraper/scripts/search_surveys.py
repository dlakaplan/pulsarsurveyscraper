#!/usr/bin/env python
import argparse
import json
import os
import re
import sys


import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from requests import Request

import pulsarsurveyscraper
import pulsarsurveyscraper.output


import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from requests import Request

import pulsarsurveyscraper
import csv
import numpy
import chime_fn
import pulsarsurveyscraper.output
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
        "--dedup",
        const=True,
        default=False,
        nargs="?",
        choices=[False, True, "hide"],
        help="Deduplicate search results?",
    )
    parser.add_argument(
        "-g",
        "--galactic",
        default=False,
        action="store_true",
        help="Search in Galactic coordinates?",
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
        "-j", "--json", default=False, action="store_true", help="Return JSON?"
    )
    parser.add_argument(
        "-v", "--verbosity", default=0, action="count", help="Increase output verbosity"
    )
    parser.add_argument(
        "--pdf",
        action="store_true",
        default=False,
        help="Save the results as pdf file?",
    )
    parser.add_argument(
        "--png",
        action="store_true",
        default=False,
        help="Save the results as png file?",
    )

    args = parser.parse_args()
    if args.verbosity == 1:
        pulsarsurveyscraper.log.setLevel("INFO")
    elif args.verbosity >= 2:
        pulsarsurveyscraper.log.setLevel("DEBUG")

    # Return to the search algorithm

    if len(args.coord) == 2:
        ra, dec = args.coord
    elif len(args.coord) == 1:
        c = args.coord[0].split()
        l = len(c)
        ra = " ".join(c[: (l // 2)])
        dec = " ".join(c[(l // 2) :])

    if not args.galactic:
        coord = pulsarsurveyscraper.parse_equcoord(ra, dec)
    else:
        coord = pulsarsurveyscraper.parse_galcoord(ra, dec)

    if coord is None:
        sys.exit(1)

    pulsar_table = pulsarsurveyscraper.PulsarTable(directory=args.dest)
    query_dict = {"type": "search",
                  "radius": args.radius}
    if not args.galactic:
        search_query_txt = "Searching {:.1f}deg around RA,Dec = {} = {}d,{}d".format(
            args.radius,
            coord.to_string("hmsdms", sep=":"),
            coord.ra.to_string(decimal=True),
            coord.dec.to_string(decimal=True, alwayssign=True),
        )
        query_dict["ra"] = coord.ra.to_string(decimal=True)
        query_dict["dec"] = coord.dec.to_string(decimal=True, alwayssign=True)

    else:
        search_query_txt = "Searching {:.1f}deg around l,b = {}d,{}d".format(
            args.radius,
            coord.l.to_string(decimal=True),
            coord.b.to_string(decimal=True, alwayssign=True),
        )
        query_dict["l"] = coord.l.to_string(decimal=True)
        query_dict["b"] = coord.b.to_string(decimal=True, alwayssign=True)


    print(search_query_txt)
    # Create the api link to be appended at the end of PDF for
    # reproducibility
    # do it without actually querying
    query_url = "https://pulsar.cgca-hub.org/api"
    response = Request("GET", query_url, params=query_dict).prepare()

    if args.dm is not None and args.dmtol is not None:
        print(
            "Also requiring DM with +/-{:.1f} of {:.1f} pc/cm**2".format(
                args.dmtol, args.dm
            )
        )
        query_dict["dm"] = args.dm
        query_dict["dmtol"] = args.dmtol

    result = pulsar_table.search(
        coord,
        radius=args.radius * u.deg,
        DM=args.dm,
        DM_tolerance=args.dmtol,
        return_json=args.json,
        deduplicate=args.dedup,
    )

    if not args.json:
        print("Found {} matches:".format(len(result)))
    else:
        print("Found {} matches:".format(result["nmatches"]))
        result = json.dumps(result)
    print(result)

    if args.png or args.pdf:
        format = "pdf" if args.pdf else "png"
        output = pulsarsurveyscraper.output.make_output(
            result, format, search_query_txt, response.url
        )
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
