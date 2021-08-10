#!/usr/bin/env python
import sys
import os
import argparse
import pulsarsurveyscraper


def main():

    all_surveys = ["all"] + list(pulsarsurveyscraper.Surveys.keys())

    pulsarsurveyscraper.log.setLevel("WARNING")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-s",
        "--survey",
        choices=all_surveys,
        nargs="*",
        help="Survey(s) to cache",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        dest="dest",
        default="",
        help="Output directory for survey caches",
    )
    parser.add_argument(
        "-v", "--verbosity", default=0, action="count", help="Increase output verbosity"
    )

    args = parser.parse_args()
    if args.verbosity == 1:
        pulsarsurveyscraper.log.setLevel("INFO")
    elif args.verbosity >= 2:
        pulsarsurveyscraper.log.setLevel("DEBUG")

    if args.survey is None or len(args.survey) == 0:
        print("Possible surveys: {}".format(", ".join(all_surveys)))
        sys.exit(0)
    if "all" in args.survey:
        surveys = pulsarsurveyscraper.Surveys.keys()
    else:
        surveys = set(args.survey)
    for survey in surveys:
        try:
            out = pulsarsurveyscraper.PulsarSurvey.read(
                survey_name=survey,
                survey_specs=pulsarsurveyscraper.Surveys[survey],
            )
            if out is None or out.data is None:
                pulsarsurveyscraper.log.error(
                    "Did not load any data for survey '{}'".format(survey)
                )
                continue
            outfile = os.path.join(args.dest, "{}.hdf5".format(survey))
            out.write(outfile, overwrite=True)
        except:
            pulsarsurveyscraper.log.error("Error loading '{}' from '{}'".format(survey,pulsarsurveyscraper.Surveys[survey]['url']))
            


if __name__ == "__main__":
    main()
