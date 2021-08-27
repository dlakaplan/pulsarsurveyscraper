#!/usr/bin/env python

#!/usr/bin/env python
import sys
import os
import argparse
import json
import re
from astropy.coordinates import SkyCoord
from astropy import units as u
import pulsarsurveyscraper
import csv
import numpy as np
import chime_fn
import logging

def main():
    pulsarsurveyscraper.log.setLevel("WARNING")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--fn", default='my_new_sources.csv', help="Clustering csv file")
    parser.add_argument("--dmtol", default=10, type=float, help="DM tolerance")
    parser.add_argument(
        "-r", "--radius", default=5, type=float, help="Search radius (degrees)"
    )
    pulsar_table = pulsarsurveyscraper.PulsarTable()

    args = parser.parse_args()
    chime_csv = chime_fn.load_new_sources(args.fn)
    unique_survey = []
    matched_survey = []
    for source in chime_csv:
        source_dict = chime_csv[source]
        ra=source_dict[0]
        dec=source_dict[1]
        dm=float(source_dict[2])
        coord = pulsarsurveyscraper.parse_equcoord(ra, dec)
        result = pulsar_table.search(
                coord,
                radius=args.radius * u.deg,
                DM=dm,
                DM_tolerance=args.dmtol,
                return_json=True,
                deduplicate=False,
        )
        if result["nmatches"]==0:
            unique_survey.append([source,float(ra),float(dec),dm])
        if result["nmatches"]>0:
            result['cluster']=[source,float(ra),float(dec),dm]
            matched_survey.append(result)
            print(result)
            print([source,float(ra),float(dec),dm])
    #now gotta check if we clash with any of the sources already scheduled by CHIME
    chime_cands = chime_fn.LoadChimeCands()
    for chime_cand in chime_cands:
        cand_coord = chime_cands[chime_cand]
        cand_ra = cand_coord[0]
        cand_dec = cand_coord[1]
        cand_dm = cand_coord[2]
        #calculate euclidean seperation
        for i,un in enumerate(unique_survey):
            euclidean_dist = np.sqrt((un[1]-cand_ra)**2+(un[2]-cand_dec)**2)
            dm_dist = np.abs(cand_dm-un[3])
            if (euclidean_dist<args.radius)&(dm_dist<args.dmtol):
                unique_survey.remove(un)
    np.save('new_sources',unique_survey)
    np.save('matched_surveys',matched_survey)

if __name__ == "__main__":
    main()
