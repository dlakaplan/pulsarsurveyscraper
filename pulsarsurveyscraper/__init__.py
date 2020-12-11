from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table, Column, vstack
from astropy.io import ascii
import numpy as np
import urllib
import json
import requests
import re
from bs4 import BeautifulSoup
from astropy import log
import argparse
import sys
import os
import time

"""
Definitions for each survey
Includes:
* the main URL
* the type (HTML, JSON, ATNF)
* period_units (ms or s)
* start_row
* pulsar_column, period_column, DM_column, ra_column, dec_column (the last two optional)
"""
# update this as needed
ATNF_version = "1.64"
Surveys = {
    "AO327": {
        "url": "http://www.naic.edu/~deneva/drift-search/index.html",
        "type": "HTML",
        "pulsar_column": 1,
        "period_column": 2,
        "DM_column": 3,
        "start_row": 1,
        "period_units": "ms",
    },
    "GBNCC": {
        "url": "http://astro.phys.wvu.edu/GBNCC/",
        "type": "HTML",
        "pulsar_column": 1,
        "period_column": 2,
        "DM_column": 3,
        "start_row": 1,
        "period_units": "ms",
    },
    "GBT350": {
        "url": "http://astro.phys.wvu.edu/GBTdrift350/",
        "type": "HTML",
        "pulsar_column": 0,
        "period_column": 1,
        "DM_column": 2,
        "start_row": 1,
        "period_units": "ms",
    },
    "PALFA": {
        "url": "http://www2.naic.edu/~palfa/newpulsars/index.html",
        "type": "HTML",
        "pulsar_column": 1,
        "period_column": 2,
        "DM_column": 3,
        "start_row": 1,
        "period_units": "ms",
    },
    "DMB": {
        "url": "http://astro.phys.wvu.edu/dmb",
        "type": "HTML",
        "pulsar_column": 0,
        "period_column": 1,
        "DM_column": 2,
        "start_row": 0,
        "period_units": "ms",
    },
    "SUPERB": {
        "url": "https://sites.google.com/site/publicsuperb/discoveries",
        "type": "HTML",
        "pulsar_column": 0,
        "period_column": 3,
        "DM_column": 4,
        "start_row": 1,
        "period_units": "ms",
        "ra_column": 1,
        "dec_column": 2,
    },
    "HTRU-S Low-latitude": {
        "url": "https://sites.google.com/site/htrusouthdeep/home/discoveries",
        "type": "HTML",
        "pulsar_column": 1,
        "period_column": 8,
        "DM_column": 9,
        "start_row": 4,
        "period_units": "ms",
        "ra_column": 6,
        "dec_column": 7,
    },
    "LOTAAS": {
        "url": "http://old.astron.nl/lotaas/index.php?sort=1",
        "type": "HTML",
        "pulsar_column": 1,
        "period_column": 2,
        "DM_column": 3,
        "start_row": 1,
        "period_units": "ms",
    },
    "RRATalog": {
        "url": "http://astro.phys.wvu.edu/rratalog",
        "type": "HTML",
        "pulsar_column": 0,
        "period_column": 1,
        "DM_column": 3,
        "start_row": 2,
        "period_units": "ms",
    },
    "CHIME": {
        "url": "http://catalog.chime-frb.ca/galactic",
        "type": "JSON",
        "period_units": "s",
    },
    "FAST": {
        "url": "http://crafts.bao.ac.cn/pulsar/",
        "type": "HTML",
        "pulsar_column": 1,
        "period_column": 4,
        "DM_column": 5,
        "start_row": 1,
        "period_units": "ms",
        "ra_column": 2,
        "dec_column": 3,
    },
    "FAST-GPPS": {
        "url": "http://zmtt.bao.ac.cn/GPPS",
        "type": "HTML",
        "pulsar_column": 0,
        "period_column": 2,
        "DM_column": 3,
        "start_row": 1,
        "period_units": "s",
        "ra_column": 4,
        "dec_column": 5,
    },
    "MWA": {
        "url": "https://wiki.mwatelescope.org/display/MP/SMART+survey+candidates",
        "type": "HTML",
        "pulsar_column": 0,
        "period_column": 1,
        "DM_column": 2,
        "start_row": 1,
        "period_units": "ms",
    },
    "ATNF": {
        "url": "https://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version={}&Name=Name&RaJ=RaJ&DecJ=DecJ&P0=P0&DM=DM&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Short+without+errors&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=51&table_bottom.y=23".format(
            ATNF_version
        ),
        "type": "ATNF",
        "period_units": "s",
    },
}


def name_to_position(name):
    name = re.sub(r"[^\d\.\+-]", "", name)
    if "-" in name:
        try:
            ra, dec = name.split("-")
        except ValueError:
            log.error("Error converting pulsar name '{}' to RA/Dec".format(name))
            return None
        sign = "-"
    else:
        ra, dec = name.split("+")
        sign = "+"
    match = re.match(
        r"J?(?P<hour>\d{2})(?P<minute>\d{2,4})(?P<decimal>\.?)(?P<frac>\d*)", ra
    )
    if match:
        if len(match.group("minute")) == 2:
            # HHMM
            ra_hms = "{}:{}{}{}".format(
                match.group("hour"),
                match.group("minute"),
                match.group("decimal"),
                match.group("frac"),
            )
        elif len(match.group("minute")) == 4:
            # HHMMSS
            ra_hms = "{}:{}:{}{}{}".format(
                match.group("hour"),
                match.group("minute")[:2],
                match.group("minute")[2:4],
                match.group("decimal"),
                match.group("frac"),
            )
        else:
            log.error("Cannot parse RA string '{}' from source '{}'".format(ra, name))
            return None
    else:
        log.error("Cannot parse RA string '{}' from source '{}'".format(ra, name))
        return None
    match = re.match(
        r"(?P<degree>\d{2})(?P<minute>\d{0,4})(?P<decimal>\.?)(?P<frac>\d*)", dec
    )
    if match:
        if len(match.group("minute")) == 0:
            # DD.D
            dec_dms = "{}{}{}".format(
                match.group("degree"), match.group("decimal"), match.group("frac"),
            )

        elif len(match.group("minute")) == 2:
            # DDMM
            dec_dms = "{}:{}{}{}".format(
                match.group("degree"),
                match.group("minute"),
                match.group("decimal"),
                match.group("frac"),
            )
        elif len(match.group("minute")) == 4:
            # DDMMSS
            dec_dms = "{}:{}:{}{}{}".format(
                match.group("degree"),
                match.group("minute")[:2],
                match.group("minute")[2:4],
                match.group("decimal"),
                match.group("frac"),
            )
        else:
            log.error("Cannot parse Dec string '{}' from source '{}'".format(dec, name))
            return None
    else:
        log.error("Cannot parse Dec string '{}' from source '{}'".format(dec, name))
        return None

    c = SkyCoord(ra_hms, sign + dec_dms, unit=("hour", "deg"))
    return c


class PulsarSurvey:
    subclasses = {}

    def __init__(
        self, survey_name=None, survey_specs=None,
    ):
        self.survey_name = survey_name
        self.survey_url = survey_specs["url"]
        self.update = None
        self.data = None
        pulsar_column = None
        period_column = None
        DM_column = None
        period_units = None
        start_row = None
        ra_column = None
        dec_column = None
        if "pulsar_column" in survey_specs:
            pulsar_column = survey_specs["pulsar_column"]
        if "period_column" in survey_specs:
            period_column = survey_specs["period_column"]
        if "DM_column" in survey_specs:
            DM_column = survey_specs["DM_column"]
        if "period_units" in survey_specs:
            period_units = survey_specs["period_units"]
        if "ra_column" in survey_specs:
            ra_column = survey_specs["ra_column"]
        if "dec_column" in survey_specs:
            dec_column = survey_specs["dec_column"]
        if "start_row" in survey_specs:
            start_row = survey_specs["start_row"]

    @classmethod
    def register(cls, survey_type):
        def decorator(subclass):
            cls.subclasses[survey_type] = subclass
            return subclass

        return decorator

    @classmethod
    def read(cls, survey_name=None, survey_specs=None):
        if survey_specs["type"] not in cls.subclasses:
            log.error("Bad survey type {}".format(survey_specs["type"]))
            return None

        return cls.subclasses[survey_specs["type"]](
            survey_name=survey_name, survey_specs=survey_specs
        )

    def write(self, filename, overwrite=False):
        self.data.write(filename, serialize_meta=True, path="data", overwrite=overwrite)
        log.info(
            "Wrote data for survey '{}' to '{}'".format(self.survey_name, filename)
        )


@PulsarSurvey.register("HTML")
class HTMLPulsarSurvey(PulsarSurvey):
    def __init__(
        self, survey_name=None, survey_specs=None,
    ):
        self.survey_name = survey_name
        self.survey_url = survey_specs["url"]
        self.update = None
        self.data = None
        pulsar_column = None
        period_column = None
        DM_column = None
        period_units = None
        start_row = None
        ra_column = None
        dec_column = None
        if "pulsar_column" in survey_specs:
            pulsar_column = survey_specs["pulsar_column"]
        if "period_column" in survey_specs:
            period_column = survey_specs["period_column"]
        if "DM_column" in survey_specs:
            DM_column = survey_specs["DM_column"]
        if "period_units" in survey_specs:
            period_units = survey_specs["period_units"]
        if "ra_column" in survey_specs:
            ra_column = survey_specs["ra_column"]
        if "dec_column" in survey_specs:
            dec_column = survey_specs["dec_column"]
        if "start_row" in survey_specs:
            start_row = survey_specs["start_row"]

        start_time = time.time()
        # parse as a HTML table
        try:
            self.page = requests.get(self.survey_url)
        except requests.exceptions.ConnectionError:
            log.error("Unable to read URL '{}'".format(self.survey_url))
            return
        self.update = Time.now()
        self.soup = BeautifulSoup(self.page.content, "html.parser")
        tables = self.soup.find_all(name="table")
        if self.survey_name == "SUPERB":
            # so far this works for SUPERB
            self.raw_table = tables[1].find(name="tr")
        elif self.survey_name == "HTRU-S Low-latitude":
            self.raw_table = tables[1]
        elif self.survey_name == "DMB" or self.survey_name == "GBT350":
            # this is very hacky
            # but the HTML seems to be missing /tr tags which breaks the parsing
            s = str(tables[0].findChildren("tr")[0])
            sout = s.replace("<tr>", "</tr><tr>")[5:].replace("<td>", "</td><td>")
            soup2 = BeautifulSoup(sout, "html.parser")
            self.raw_table = soup2
        else:
            self.raw_table = tables[0]

        self.rows = self.raw_table.find_all(name="tr")
        pulsar = []
        period = []
        DM = []
        RA = []
        Dec = []
        for row in self.rows[start_row:]:
            cols = row.find_all(name="td")
            if (
                (len(cols) < 3)
                or ("pulsar" in cols[pulsar_column].text.lower())
                or ("name" in cols[pulsar_column].text.lower())
            ):
                continue
            name = cols[pulsar_column].text
            # replace some dashes with minus signs
            name = name.replace(chr(8211), "-")
            name = re.sub(r"[^J\d\+-\.A-Za-g]", "", name)
            if name.startswith("FRB") or len(name) == 0:
                continue
            pulsar.append(name.strip())
            P = cols[period_column].text
            # special cases and unit conversion
            P = re.sub(r"[^\d\.]", "", P)
            if period_units == "ms":
                try:
                    period.append(float(P))
                except ValueError:
                    period.append(np.nan)
            elif period_units == "s":
                period.append(float(P) * 1000)
            try:
                dm = re.sub(r"[^\d\.]", "", cols[DM_column].text)
                DM.append(float(dm))
            except ValueError:
                log.error(
                    "Error parsing DM value of '{}' for pulsar '{}'".format(
                        cols[DM_column].text, pulsar[-1]
                    )
                )
                return
            if ra_column is None or dec_column is None:
                try:
                    coord = name_to_position(pulsar[-1])
                except:
                    log.warning(
                        "Unable to parse pulsar '{}' to determine coordiates; assuming (0,0)".format(
                            pulsar[-1]
                        )
                    )
                    coord = SkyCoord(0 * u.deg, 0 * u.deg)
            else:
                ra_text = re.sub(r"[^\d:\.]", "", cols[ra_column].text)
                dec_text = re.sub(r"[^\d:\.\+-]", "", cols[dec_column].text)
                if len(ra_text) == 0 or len(dec_text) == 0:
                    try:
                        coord = name_to_position(pulsar[-1])
                    except:
                        log.warning(
                            "No RA/Dec available and unable to parse pulsar '{}' to determine coordiates; assuming (0,0)".format(
                                pulsar[-1]
                            )
                        )
                        coord = SkyCoord(0 * u.deg, 0 * u.deg)
                else:
                    try:
                        coord = SkyCoord(ra_text, dec_text, unit=("hour", "deg"),)
                    except ValueError:
                        log.error(
                            "Error parsing position values of '{},{}' for pulsar '{}'".format(
                                cols[ra_column].text, cols[dec_column].text, pulsar[-1],
                            )
                        )
                        return
            RA.append(coord.ra.deg)
            Dec.append(coord.dec.deg)
        self.data = Table(
            [
                Column(pulsar, name="PSR"),
                Column(RA, name="RA", unit=u.deg, format="%.6f"),
                Column(Dec, name="Dec", unit=u.deg, format="%.6f"),
                Column(period, name="P", unit=u.ms, format="%.2f"),
                Column(DM, name="DM", unit=u.pc / u.cm ** 3, format="%.2f"),
            ]
        )
        end_time = time.time()
        log.info(
            "Read data for {} pulsars for survey '{}' in {:.2f}s at {}".format(
                len(self.data), self.survey_name, end_time - start_time, self.update.iso
            )
        )
        self.data.meta["url"] = self.survey_url
        self.data.meta["survey"] = self.survey_name
        self.data.meta["date"] = self.update


@PulsarSurvey.register("ATNF")
class ATNFPulsarSurvey(PulsarSurvey):
    def __init__(
        self, survey_name=None, survey_specs=None,
    ):
        self.survey_name = survey_name
        self.survey_url = survey_specs["url"]
        self.update = None
        self.data = None
        pulsar_column = None
        period_column = None
        DM_column = None
        period_units = None
        start_row = None
        ra_column = None
        dec_column = None
        if "pulsar_column" in survey_specs:
            pulsar_column = survey_specs["pulsar_column"]
        if "period_column" in survey_specs:
            period_column = survey_specs["period_column"]
        if "DM_column" in survey_specs:
            DM_column = survey_specs["DM_column"]
        if "period_units" in survey_specs:
            period_units = survey_specs["period_units"]
        if "ra_column" in survey_specs:
            ra_column = survey_specs["ra_column"]
        if "dec_column" in survey_specs:
            dec_column = survey_specs["dec_column"]
        if "start_row" in survey_specs:
            start_row = survey_specs["start_row"]

        if self.survey_name == "ATNF":
            start_time = time.time()
            try:
                self.page = requests.get(self.survey_url)
            except requests.exceptions.ConnectionError:
                log.error("Unable to read URL '{}'".format(self.survey_url))
                return
            self.update = Time.now()
            self.soup = BeautifulSoup(self.page.content, "html.parser")
            self.raw_table = self.soup.find(name="pre").text
            data = ascii.read(
                self.raw_table.replace("#", "N"),
                comment="----",
                fill_values=("*", -99),
                data_start=2,
            )
            coord = SkyCoord(data["RAJ"], data["DECJ"], unit=("hour", "deg"))
            data["P0"].unit = u.s
            data["DM"].unit = u.pc / u.cm ** 3
            data["NAME"].name = "PSR"
            self.data = Table(
                (
                    data["PSR"],
                    Column(coord.ra.deg * u.deg, name="RA", unit=u.deg, format="%.6f"),
                    Column(
                        coord.dec.deg * u.deg, name="Dec", unit=u.deg, format="%.6f"
                    ),
                    Column(data["P0"].to(u.ms), name="P", format="%.2f"),
                    data["DM"],
                )
            )
        end_time = time.time()
        log.info(
            "Read data for {} pulsars for survey '{}' in {:.2f}s at {}".format(
                len(self.data), self.survey_name, end_time - start_time, self.update.iso
            )
        )
        self.data.meta["url"] = self.survey_url
        self.data.meta["survey"] = self.survey_name
        self.data.meta["date"] = self.update


@PulsarSurvey.register("JSON")
class JSONPulsarSurvey(PulsarSurvey):
    def __init__(
        self, survey_name=None, survey_specs=None,
    ):
        self.survey_name = survey_name
        self.survey_url = survey_specs["url"]
        self.update = None
        self.data = None
        pulsar_column = None
        period_column = None
        DM_column = None
        period_units = None
        start_row = None
        ra_column = None
        dec_column = None
        if "pulsar_column" in survey_specs:
            pulsar_column = survey_specs["pulsar_column"]
        if "period_column" in survey_specs:
            period_column = survey_specs["period_column"]
        if "DM_column" in survey_specs:
            DM_column = survey_specs["DM_column"]
        if "period_units" in survey_specs:
            period_units = survey_specs["period_units"]
        if "ra_column" in survey_specs:
            ra_column = survey_specs["ra_column"]
        if "dec_column" in survey_specs:
            dec_column = survey_specs["dec_column"]
        if "start_row" in survey_specs:
            start_row = survey_specs["start_row"]

        # read as JSON
        # with urllib.request.urlopen(self.survey_url) as url:
        #    self.raw_table = json.loads(url.read().decode())
        start_time = time.time()
        req = urllib.request.Request(self.survey_url)
        try:
            response = urllib.request.urlopen(req)
        except urllib.error.URLError as e:
            if hasattr(e, "reason"):
                log.error(
                    "Failed to reach server '{}': {}".format(self.survey_url, e.reason)
                )
            elif hasattr(e, "code"):
                log.error(
                    "Server '{}' could not fulfill request: code {}".format(
                        self.survey_url, e.code
                    )
                )
            return
        self.raw_table = json.loads(response.read().decode())
        self.update = Time.now()
        pulsar = []
        period = []
        DM = []
        RA = []
        Dec = []
        for key in self.raw_table.keys():
            pulsar.append(key)
            coord = SkyCoord(
                self.raw_table[key]["ra"]["value"],
                self.raw_table[key]["dec"]["value"],
                unit=("hour", "deg"),
            )
            RA.append(coord.ra.deg)
            Dec.append(coord.dec.deg)
            if period_units == "ms":
                try:
                    period.append(float(self.raw_table[key]["period"]["value"]))
                except TypeError:
                    period.append(np.nan)
            elif period_units == "s":
                try:
                    period.append(float(self.raw_table[key]["period"]["value"]) * 1000)
                except TypeError:
                    period.append(np.nan)

            DM.append(self.raw_table[key]["dm"]["value"])

        self.data = Table(
            [
                Column(pulsar, name="PSR"),
                Column(RA, name="RA", unit=u.deg, format="%.6f"),
                Column(Dec, name="Dec", unit=u.deg, format="%.6f"),
                Column(period, name="P", unit=u.ms, format="%.2f"),
                Column(DM, name="DM", unit=u.pc / u.cm ** 3, format="%.2f"),
            ]
        )
        end_time = time.time()
        log.info(
            "Read data for {} pulsars for survey '{}' in {:.2f}s at {}".format(
                len(self.data), self.survey_name, end_time - start_time, self.update.iso
            )
        )
        self.data.meta["url"] = self.survey_url
        self.data.meta["survey"] = self.survey_name
        self.data.meta["date"] = self.update


class PulsarTable:
    def __init__(self, directory=None):
        if directory is None:
            directory = os.path.curdir
        self.directory = os.path.abspath(directory)
        data = []
        for survey in Surveys:
            surveyfile = os.path.join(self.directory, "{}.hdf5".format(survey))
            if os.path.exists(surveyfile):
                log.info("Reading survey '{}' from {}".format(survey, surveyfile))
                data.append(Table.read(surveyfile, path="data"))
            else:
                log.info("Cannot find cached survey '{}': loading...".format(survey))
                out = pulsarsurveyscraper.PulsarSurvey.read(
                    survey_name=survey,
                    survey_specs=pulsarsurveyscraper.Surveys[survey],
                )
                data.append(out.data)
            data[-1].add_column(
                Column(np.array([survey] * len(data[-1])), name="survey")
            )
            data[-1].add_column(
                Column(
                    np.array([data[-1].meta["date"]] * len(data[-1])),
                    name="retrieval date",
                )
            )
            data[-1].meta = {}
        self.data = vstack(data)
        self.coord = SkyCoord(self.data["RA"], self.data["Dec"])

    def search(self, coord, radius=1 * u.deg, DM=None, DM_tolerance=10):
        """
        out = search(coord, radius=1 * u.deg, DM=None, DM_tolerance=10)

        returns an astropy.table object with the sources that match the given criteria
        """
        distance = self.coord.separation(coord)
        i = np.argsort(distance)
        distance = distance[i]
        data = self.data
        good = distance < radius
        if DM is not None and DM_tolerance is not None:
            good = good & (np.fabs(data["DM"].data - DM) < DM_tolerance)
        good = good
        output = data[good]
        output.add_column(Column(distance[good], name="Distance", format=".4f"))
        return output
