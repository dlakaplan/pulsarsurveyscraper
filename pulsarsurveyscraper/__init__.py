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

# import from another file
from .surveys import Surveys

def name_to_position(name: str) -> SkyCoord:
    """
    parses a pulsar name like J1234+5656 and returns an astropy SkyCoord object
    returns None if parsing fails

    Args:
        name (str): name to parse

    Returns:
        SkyCoord: the coordinates corresponding to the name (or None)

    formats:
    J1234+56
    J1234.5+56
    J123456+56
    J123456.7+56
    J1234+12.3
    J1234+1234
    J1234+1234.5
    J1234+123456
    J1234+123456.7

    """
    # remove any characters that are not a digit, decimal, or sign
    name = re.sub(r"[^\d\.\+-]", "", name)
    if "-" in name:
        try:
            ra, dec = name.split("-")
        except ValueError as e:
            log.error("Error converting pulsar name '{}' to RA/Dec: {}".format(name, e))
            return None
        sign = "-"
    else:
        try:
            ra, dec = name.split("+")
        except ValueError as e:
            log.error("Error converting pulsar name '{}' to RA/Dec: {}".format(name, e))
            return None
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
                match.group("degree"),
                match.group("decimal"),
                match.group("frac"),
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

    try:
        c = SkyCoord(ra_hms, sign + dec_dms, unit=("hour", "deg"))
    except ValueError as e:
        log.error("Cannot parse RA/Dec {},{}: {}".format(ra_hms, sign + dec_dms, e))
        return None
    return c


class PulsarSurvey:
    """
    PulsarSurvey()

    top-level class for pulsar survey scraping
    defines class-level read() function which goes to appropriate subclass
    also has write() function

    survey_specs contains:
    * the main URL
    * the type (HTML, JSON, ATNF)
    * period_units (ms or s)
    * start_row (optional for JSON/ATNF)
    * pulsar_column, period_column, DM_column, ra_column, dec_column (the last two optional for all)

    general logic is to read a survey from its website
    construct astropy Table with the contents

    columns:
    PSR: name of pulsar
    P: period of pulsar (ms)
    DM: DM of pulsar
    RA: RA of pulsar (deg)
    Dec: Dec of pulsar (deg)

    also adds metadata about the survey name, site URL, and date retrieved

    Args:
       survey_name (str): name of survey
       survey_specs (dict): dictionary containing other needed info (such as "url")
    """

    subclasses = {}

    def __init__(
        self,
        survey_name: str = None,
        survey_specs: dict = None,
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
    def read(cls, survey_name: str = None, survey_specs: dict = None):
        """
        Args:
            survey_name (str): name of survey
            survey_specs (dict): dictionary containing other needed info (such as "url")

        Returns:
            PulsarSurvey: instance of the appropriate subclass

        """
        if survey_specs["type"] not in cls.subclasses:
            log.error("Bad survey type {}".format(survey_specs["type"]))
            return None

        return cls.subclasses[survey_specs["type"]](
            survey_name=survey_name, survey_specs=survey_specs
        )

    def write(self, filename: str, overwrite: bool = False):
        """
        Writes data to HDF5 file

        Args:
            filename (str): name of file for output
            overwrite (bool, optional): whether to overwrite existing output files.  Defaults to False

        """
        self.data.write(filename, serialize_meta=True, path="data", overwrite=overwrite)
        log.info(
            "Wrote data for survey '{}' to '{}'".format(self.survey_name, filename)
        )


@PulsarSurvey.register("HTML")
class HTMLPulsarSurvey(PulsarSurvey):
    """
    HTMLPulsarSurvey(PulsarSurvey)

    subclass appropriate for pages contained in HTML <table>s
    some special logic built-in for handling google tables (SUPERB)
    and poorly-formatted HTML (DMB, GBT350)

    """

    def __init__(
        self,
        survey_name: str = None,
        survey_specs: dict = None,
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
            # iterate over each row in the table
            # each row represents a pulsar (usually)
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
            except ValueError as e:
                log.error(
                    "Error parsing DM value of '{}' for pulsar '{}': {}".format(
                        cols[DM_column].text, pulsar[-1], e
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
                        coord = SkyCoord(
                            ra_text,
                            dec_text,
                            unit=("hour", "deg"),
                        )
                    except ValueError as e:
                        log.error(
                            "Error parsing position values of '{},{}' for pulsar '{}': {}".format(
                                cols[ra_column].text,
                                cols[dec_column].text,
                                pulsar[-1],
                                e,
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
    """
    ATNFPulsarSurvey(PulsarSurvey)

    subclass appropriate for ATNF pulsar database
    which returns plain-text table (no HTML)
    """

    def __init__(
        self,
        survey_name: str = None,
        survey_specs: dict = None,
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
    """
    JSONPulsarSurvey(PulsarSurvey)


    subclass for surveys that return JSON objects
    currently JSON keys only defined for CHIME
    """

    def __init__(
        self,
        survey_name: str = None,
        survey_specs: dict = None,
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
    """
    PulsarTable

    construct a table containing the results of the individual surveys, either from cached copies or by direct queries

    adds columns to table with:
    survey: survey name
    retrieval date: date survey content was retrieved

    has method search() which allows cone searches (with or without DM constraints)

    Args:
        directory (str, optional): directory to search for cached files.  Defaults to current directory

    """

    def __init__(self, directory: str = None):
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
                log.info(
                    "Cannot find cached survey '{}' at {}: loading...".format(
                        survey, surveyfile
                    )
                )
                out = PulsarSurvey.read(
                    survey_name=survey,
                    survey_specs=Surveys[survey],
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

    def search(
        self,
        coord: SkyCoord,
        radius: u.quantity.Quantity = 1 * u.deg,
        DM: float = None,
        DM_tolerance: float = 10,
        return_json: bool = False,
    ):
        """
        returns an astropy.table object with the sources that match the given criteria
        or a JSON object if return_json=True

        Args:
            coord (SkyCoord): position to search around
            radius (u.quantity.Quantity, optional): radius to search
            DM (float, optional): central DM constraint (if None, do no DM cut)
            DM_tolerance (float, optional): DM tolerance for constraint
            return_json (bool, optional): return JSON instead of Table

        Returns:
            Table or dict
        """
        distance = self.coord.separation(coord)
        i = np.argsort(distance)
        distance = distance[i]
        data = self.data[i]
        good = distance < radius
        if DM is not None and DM_tolerance is not None:
            good = good & (np.fabs(data["DM"].data - DM) < DM_tolerance)
        good = good
        output = data[good]
        output.add_column(Column(distance[good], name="Distance", format=".4f"))
        if not return_json:
            return output
        output_dict = {}
        output_dict["searchra"] = {
            "display_name": "Search RA (deg)",
            "value": coord.ra.deg,
        }
        output_dict["searchdec"] = {
            "display_name": "Search Dec (deg)",
            "value": coord.dec.deg,
        }
        output_dict["searchcoord"] = {
            "display_name": "Search Coord",
            "value": coord.to_string("hmsdms", sep=":"),
        }
        output_dict["searchrad"] = {
            "display_name": "Search Radius (deg)",
            "value": radius.value,
        }
        if DM is not None and DM_tolerance is not None:
            output_dict["searchdm"] = {"display_name": "Search DM", "value": DM}
            output_dict["searchdmtolerance"] = {
                "display_name": "DM Tolerance",
                "value": DM_tolerance,
            }

        output_dict["nmatches"] = len(output)
        for row in output:
            key = row["PSR"]
            output_dict[key] = {
                "ra": {"display_name": "RA (deg)", "value": row["RA"]},
                "dec": {"display_name": "Dec (deg)", "value": row["Dec"]},
                "period": {"display_name": "Spin Period (ms)", "value": row["P"]},
                "dm": {"display_name": "DM (pc/cc)", "value": row["DM"]},
                "survey": {"display_name": "Survey", "value": row["survey"]},
                "url": {
                    "display_name": "Survey URL",
                    "value": Surveys[row["survey"]]["url"],
                },
                "distance": {
                    "display_name": "Distance (deg)",
                    "value": row["Distance"],
                },
                "date": {
                    "display_name": "Retrieval Date",
                    "value": row["retrieval date"].iso,
                },
            }
        return json.dumps(output_dict)
