import flask
from flask import request, jsonify
import sys
import os
import argparse
import re
from astropy.coordinates import SkyCoord
from astropy import units as u
import pulsarsurveyscraper
import pandas

app = flask.Flask(__name__)
app.config["DEBUG"] = True


@app.route("/", methods=["GET"])
def home():
    return """<h1>Pulsar Survey Scraper</h1>
<p>A prototype API for searching pulsar surveys.</p>"""


@app.route("/api/v1/surveys/all", methods=["GET"])
def api_all():
    """
    return list of all surveys
    """
    surveys = {"surveys": list(pulsarsurveyscraper.Surveys.keys())}
    return jsonify(surveys)


@app.route("/api/v1/surveys/search", methods=["GET"])
def api_search(radius=5, dm=None, dmtol=10):

    if "ra" in request.args:
        ra = float(request.args["ra"])
    else:
        return "Error: no RA specified"
    if "dec" in request.args:
        dec = float(request.args["dec"])
    else:
        return "Error: no Dec specified"
    if "radius" in request.args:
        radius = float(request.args["radius"])
    if "dm" in request.args:
        dm = float(request.args["dm"])
    if "dmtol" in request.args:
        dmtol = float(request.args["dmtol"])

    coord = SkyCoord(ra * u.deg, dec * u.deg)
    result = pulsar_table.search(
        coord, radius=radius * u.deg, DM=dm, DM_tolerance=dmtol, return_json=True
    )

    return result


@app.route("/api/v1/surveys/display", methods=["GET"])
def api_display(radius=5, dm=None, dmtol=10):

    if "ra" in request.args:
        ra = float(request.args["ra"])
    else:
        return "Error: no RA specified"
    if "dec" in request.args:
        dec = float(request.args["dec"])
    else:
        return "Error: no Dec specified"
    if "radius" in request.args:
        radius = float(request.args["radius"])
    if "dm" in request.args:
        dm = float(request.args["dm"])
    if "dmtol" in request.args:
        dmtol = float(request.args["dmtol"])

    coord = SkyCoord(ra * u.deg, dec * u.deg)
    result = pulsar_table.search(
        coord, radius=radius * u.deg, DM=dm, DM_tolerance=dmtol
    )
    df = result.to_pandas()

    result_string = "<h1>Pulsar Survey Scraper</h1>\n"
    result_string += "<h3>Searching {:.1f}d around RA,Dec {} = {}d,{}d</h3>\n".format(
        radius,
        coord.to_string("hmsdms", sep=":"),
        coord.ra.to_string(decimal=True),
        coord.dec.to_string(decimal=True, alwayssign=True),
    )
    if dm is not None and dmtol is not None:
        result_string += "<h3>Also requiring DM with +/-{:.1f} of {:.1f} pc/cm**2</h3>\n".format(
            dmtol, dm
        )
    # originally this is a byte-string
    # reformat as a normal string
    df["PSR"] = df["PSR"].str.decode("utf-8")
    html_table = df.to_html(
        formatters={
            "P": lambda x: "{:<.2f}".format(x),
            "Distance": lambda x: "%.2f" % x,
        }
    )
    # reformat a bit to get links to the survey sites
    for survey in pulsarsurveyscraper.Surveys:
        html_table = html_table.replace(
            survey,
            "<a href='{}'>{}</a>".format(
                pulsarsurveyscraper.Surveys[survey]["url"], survey
            ),
        )
    result_string += html_table
    return result_string


def main():
    pulsarsurveyscraper.log.setLevel("WARNING")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
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

    # how can I avoid making this a global?
    global pulsar_table
    pulsar_table = pulsarsurveyscraper.PulsarTable(directory=args.dest)

    app.run()


if __name__ == "__main__":
    main()
