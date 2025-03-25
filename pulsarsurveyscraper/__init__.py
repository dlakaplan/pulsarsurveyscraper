import argparse
import json
import os
import re
import sys
import time
import typing
import urllib

import numpy as np
import requests
from astropy import log
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Column, Table, MaskedColumn, vstack
from astropy.time import Time
from bs4 import BeautifulSoup

# import from another file
from .surveys import Surveys


def extract_from_json(jsondict, keylist):
    """
    Extract a value from a JSON dictionary and a list of keys
    """
    if len(keylist) == 1:
        return jsondict[keylist[0]]
    return extract_from_json(jsondict[keylist[0]], keylist[1:])


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
    # also replace unicode 8722 (minus sign) with standard "-"
    name = re.sub(chr(8722), "-", name)
    name = re.sub(r"[^\d\.\+-]", "", name)
    if not ("+" in name or "-" in name):
        log.error(
            f"Error converting pulsar name '{name}' to RA/Dec: name {name} does not contain '+' or '-'"
        )
        return None
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


def parse_equcoord(ra, dec):
    """
    parse_equcoord(ra,dec)

    parse strings ra,dec and try to turn into a SkyCoord

    Args:
        ra (str): ra string
        dec (str): dec string

    Returns:
        SkyCoord or None (on failure)

    """
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
            log.error("Unable to parse input coordinates '{},{}'".format(ra, dec))
            return None
    return coord


def parse_galcoord(l, b):
    """
    parse_equcoord(l,b)

    parse strings l,b and try to turn into a SkyCoord

    Args:
        l (str): l string
        b (str): b string

    Returns:
        SkyCoord or None (on failure)

    """
    try:
        if (re.search(r"[^\d.+\-]", l) is None) and (
            re.search(r"[^\d.+\-]", b) is None
        ):
            coord = SkyCoord(l, b, unit="deg", frame="galactic")
        else:
            coord = SkyCoord(l, b, frame="galactic")
    except ValueError:
        log.error("Unable to parse input coordinates '{},{}'".format(ra, dec))
        return None
    return coord


def deduplicate_table(
    data,
    period_error=0.02,
    DM_error=0.02,
    distance=1 * u.deg,
    precedence={"ATNF": 0, "GalacticMSPs": 1},
):
    """Deduplicate the pulsar table

    Identifies potential duplicates based on position match, period match, DM match
    does not remove them but flags them in the table (adds a column named "Duplicate?")

    Parameters
    ----------
    data : astropy.table.Table
        Input pulsar table
    period_error : float
        fractional period tolerance for a match (default = 0.02)
    DM_error : float
        fractional DM tolerance for a mtch (default = 0.02)
    distance: astropy.unit.Quantity
        angular distance for a match (default = 1*u.deg)
    precedence: dict
        precedence order for surveys.  If not listed all others are equivalent

    Adds a column to self.data named "Duplicate?" with potential matches
    """

    coord = SkyCoord(data["RA"], data["Dec"])

    i1, i2, d2, _ = coord.search_around_sky(coord, seplimit=distance)
    unique_sources = []
    duplicate_sources = []
    duplications = {}
    duplicate_column = [None] * len(data)

    for j in i1:
        if j in duplicate_sources:
            # already identified as a duplicate
            continue
        possible_matches = i2[(i1 == j) & (i2 != j)]
        # make sure the potential matches weren't already identified
        possible_matches = np.setdiff1d(possible_matches, unique_sources)
        possible_matches = np.setdiff1d(possible_matches, duplicate_sources)
        if possible_matches.sum() == 0:
            unique_sources.append(j)
            continue
        period_match = (
            np.fabs((data["P"][possible_matches] - data["P"][j]) / data["P"][j])
            < period_error
        )
        DM_match = (
            np.fabs((data["DM"][possible_matches] - data["DM"][j]) / data["DM"][j])
            < DM_error
        )
        # we want ones that end is a number (not a letter)
        isnotGC = np.array(
            [re.search(r"\d$", x) is not None for x in data["PSR"][possible_matches]]
        )

        match = period_match & DM_match & isnotGC

        if match.any():
            all_matches = np.append(possible_matches[match], j)
            order = np.zeros(len(all_matches))
            for k in range(len(all_matches)):
                if data["survey"][all_matches[k]] in precedence:
                    order[k] = precedence[data["survey"][all_matches[k]]]
                else:
                    order[k] = len(precedence)
            all_matches = all_matches[np.argsort(order)]
            unique = all_matches[0]
            dups = all_matches[1:]
            unique_sources.append(unique)
            duplicate_sources += list(dups)
            duplications[unique] = list(dups)
            for k in dups:
                duplicate_column[k] = f"{data[unique]['survey']}:{data[unique]['PSR']}"
                # print(f"Found sources {list(dups)} are duplicates of {unique}")
                # print(data[all_matches])
            else:
                unique_sources.append(j)

    data.add_column(Column(duplicate_column, name="Duplicate?"))


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
    * coordinate_frame, ra_unit, dec_unit (if ra_column/dec_column supplied)
    * table_index (the table number on a page)

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

    has method load_specs() which sets default parameter values and overides if others are supplied
    """

    subclasses = {}

    def __init__(
        self,
        survey_name: str = None,
        survey_specs: dict = None,
    ):
        self.survey_name = survey_name
        self.load_specs(survey_specs)

    @classmethod
    def register(cls, survey_type):
        def decorator(subclass):
            cls.subclasses[survey_type] = subclass
            return subclass

        return decorator

    def load_specs(self, survey_specs):
        """
        Args:
            survey_specs (dict): dictionary with parameter values to override defaults

        """
        # defaults
        self.update = None
        self.data = None
        self.pulsar_column = None
        self.period_column = None
        self.DM_column = None
        self.period_units = None
        self.start_row = None
        self.ra_column = None
        self.dec_column = None
        self.coordinate_frame = "icrs"
        self.ra_unit = "hour"
        self.dec_unit = "deg"
        self.table_index = 0
        self.survey_url = survey_specs["url"]
        for k in survey_specs:
            self.__dict__[k] = survey_specs[k]

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
        self.load_specs(survey_specs)

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
            self.raw_table = tables[self.table_index]

        self.rows = self.raw_table.find_all(name="tr")
        pulsar = []
        period = []
        DM = []
        RA = []
        Dec = []
        for row in self.rows[self.start_row :]:
            # iterate over each row in the table
            # each row represents a pulsar (usually)
            cols = row.find_all(name="td")
            if (
                (len(cols) < 3)
                or ("pulsar" in cols[self.pulsar_column].text.lower())
                or ("name" in cols[self.pulsar_column].text.lower())
            ):
                continue
            name = cols[self.pulsar_column].text
            # replace some dashes with minus signs
            name = name.replace(chr(8211), "-")
            name = name.replace(chr(8722), "-")
            name = re.sub(r"[^J\d\+-\.A-Za-z]", "", name)
            if name.startswith("FRB") or len(name) == 0:
                continue
            pulsar.append(name.strip())
            P = cols[self.period_column].text
            # get rid of errors in parentheses
            P = re.sub(r"\(\S+\)", "", P)
            # special cases and unit conversion
            P = re.sub(r"[^\d\.]", "", P)
            if self.period_units == "ms":
                try:
                    period.append(float(P))
                except ValueError:
                    period.append(np.nan)
            elif self.period_units == "s":
                try:
                    period.append(float(P) * 1000)
                except ValueError:
                    period.append(np.nan)
            try:
                # get rid of errors in parentheses
                if self.survey_name == "Puschino":
                    # they have a range of DM
                    if "-" in cols[self.DM_column].text:
                        dms = list(map(float, cols[self.DM_column].text.split("-")))
                        dm = 0.5 * (dms[0] + dms[1])
                    else:
                        dm = cols[self.DM_column].text
                else:
                    dm = re.sub(r"\(\S+\)", "", cols[self.DM_column].text)
                    dm = re.sub(r"[^\d\.]", "", dm)
                DM.append(float(dm))
            except ValueError as e:
                log.error(
                    "Error parsing DM value of '{}' for pulsar '{}': {}".format(
                        cols[self.DM_column].text, pulsar[-1], e
                    )
                )
                return
            if self.ra_column is None or self.dec_column is None:
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
                if self.survey_name == "Puschino":
                    ra_text = re.sub(
                        r"s$",
                        "",
                        re.sub(r"m", ":", re.sub(r"h", ":", cols[self.ra_column].text)),
                    )
                else:
                    ra_text = re.sub(r"\(\S+\)", "", cols[self.ra_column].text)
                    ra_text = re.sub(r"[^\d:\.]", "", ra_text)
                # some of the HTML tables have a non-breaking hyphen (Unicode 8209)
                # instead of a hyphen
                # convert it
                dec_text = cols[self.dec_column].text
                if chr(8209) in dec_text:
                    dec_text = dec_text.replace(chr(8209), "-")
                if self.survey_name == "Puschino":
                    dec_text = re.sub(r"\'", ":00", re.sub(r"o", ":", dec_text))
                dec_text = re.sub(r"\(\S+\)", "", dec_text)
                dec_text = re.sub(r"[^\d:\.\+-]", "", dec_text)
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
                            frame=self.coordinate_frame,
                            unit=(self.ra_unit, self.dec_unit),
                        ).icrs
                    except ValueError as e:
                        log.error(
                            "Error parsing position values of '{},{}' for pulsar '{}': {}".format(
                                cols[self.ra_column].text,
                                cols[self.dec_column].text,
                                pulsar[-1],
                                e,
                            )
                        )
                        return
            if coord is None:
                log.warning(
                    "Unable to parse pulsar '{}'; assuming (0,0).".format(pulsar[-1])
                )
                coord = SkyCoord(0 * u.deg, 0 * u.deg)
            RA.append(coord.ra.deg)
            Dec.append(coord.dec.deg)
            if (
                len(pulsar) >= 2
                and pulsar[-1] == pulsar[-2]
                and period[-1] == period[-2]
                and DM[-1] == DM[-2]
            ):
                log.warning(
                    f"Identified apparent duplicate:\n\t{pulsar[-1]} {period[-1]} {DM[-1]}\n\t{pulsar[-2]} {period[-2]} {DM[-2]}\nDeleting..."
                )
                # it's a duplicate
                del pulsar[-1]
                del period[-1]
                del DM[-1]
                del RA[-1]
                del Dec[-1]

        self.data = Table(
            [
                Column(pulsar, name="PSR"),
                Column(RA, name="RA", unit=u.deg, format="%.6f"),
                Column(Dec, name="Dec", unit=u.deg, format="%.6f"),
                Column(period, name="P", unit=u.ms, format="%.2f"),
                Column(DM, name="DM", unit=u.pc / u.cm**3, format="%.2f"),
            ]
        )
        end_time = time.time()
        log.info(
            "Read data for {} pulsars for survey '{}' in {:.2f}s at {}".format(
                len(self.data),
                self.survey_name,
                end_time - start_time,
                self.update.to_value("iso", subfmt="date_hm"),
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
    which returns plain-text table (no HTML) but in a HTML table
    """

    def __init__(
        self,
        survey_name: str = None,
        survey_specs: dict = None,
    ):
        self.survey_name = survey_name
        self.load_specs(survey_specs)

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
            data["DM"].unit = u.pc / u.cm**3
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
                len(self.data),
                self.survey_name,
                end_time - start_time,
                self.update.to_value("iso", subfmt="date_hm"),
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
        self.load_specs(survey_specs)

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
            if self.ra_key is not None and self.dec_key is not None:
                coord = SkyCoord(
                    extract_from_json(self.raw_table[key], self.ra_key),
                    extract_from_json(self.raw_table[key], self.dec_key),
                    unit=("hour", "deg"),
                )
            else:
                try:
                    coord = name_to_position(pulsar[-1])
                except:
                    log.warning(
                        "No RA/Dec available and unable to parse pulsar '{}' to determine coordiates; assuming (0,0)".format(
                            pulsar[-1]
                        )
                    )
                    coord = SkyCoord(0 * u.deg, 0 * u.deg)
            RA.append(coord.ra.deg)
            Dec.append(coord.dec.deg)
            if self.period_units == "ms":
                try:
                    period.append(
                        float(extract_from_json(self.raw_table[key], self.period_key))
                    )
                except TypeError:
                    period.append(np.nan)
            elif self.period_units == "s":
                try:
                    period.append(
                        float(extract_from_json(self.raw_table[key], self.period_key))
                        * 1000
                    )
                except TypeError:
                    period.append(np.nan)

            DM.append(extract_from_json(self.raw_table[key], self.dm_key))

        self.data = Table(
            [
                Column(pulsar, name="PSR"),
                Column(RA, name="RA", unit=u.deg, format="%.6f"),
                Column(Dec, name="Dec", unit=u.deg, format="%.6f"),
                Column(period, name="P", unit=u.ms, format="%.2f"),
                Column(DM, name="DM", unit=u.pc / u.cm**3, format="%.2f"),
            ]
        )
        end_time = time.time()
        log.info(
            "Read data for {} pulsars for survey '{}' in {:.2f}s at {}".format(
                len(self.data),
                self.survey_name,
                end_time - start_time,
                self.update.to_value("iso", subfmt="date_hm"),
            )
        )
        self.data.meta["url"] = self.survey_url
        self.data.meta["survey"] = self.survey_name
        self.data.meta["date"] = self.update


@PulsarSurvey.register("ASCII")
class ASCIIPulsarSurvey(PulsarSurvey):
    """
    ASCIIPulsarSurvey(PulsarSurvey)

    subclass appropriate for ASCII pulsar database
    which returns plain-text table (no HTML)
    """

    def __init__(
        self,
        survey_name: str = None,
        survey_specs: dict = None,
    ):
        self.survey_name = survey_name
        self.load_specs(survey_specs)

        start_time = time.time()
        try:
            self.page = requests.get(self.survey_url)
        except requests.exceptions.ConnectionError:
            log.error("Unable to read URL '{}'".format(self.survey_url))
            return
        self.update = Time.now()
        self.raw_table = self.page.content.decode("utf-8")
        data = ascii.read(
            self.raw_table,
            format="no_header",
            fill_values=(("*", -99), ("N/A", -99), ("unk", -99)),
        )
        coord = SkyCoord(
            data.columns[self.ra_column],
            data.columns[self.dec_column],
            unit=(self.ra_unit, self.dec_unit),
            frame=self.coordinate_frame,
        ).icrs
        data.columns[self.period_column].name = "P"
        if self.period_units == "s":
            data["P"].unit = u.s
        else:
            data["P"].unit = u.ms
        data.columns[self.DM_column].name = "DM"
        data["DM"].unit = u.pc / u.cm**3
        data.columns[self.pulsar_column].name = "PSR"
        self.data = Table(
            (
                data["PSR"],
                Column(coord.ra.deg * u.deg, name="RA", unit=u.deg, format="%.6f"),
                Column(coord.dec.deg * u.deg, name="Dec", unit=u.deg, format="%.6f"),
                Column(data["P"].to(u.ms), name="P", format="%.2f"),
                data["DM"],
            )
        )
        end_time = time.time()
        log.info(
            "Read data for {} pulsars for survey '{}' in {:.2f}s at {}".format(
                len(self.data),
                self.survey_name,
                end_time - start_time,
                self.update.to_value("iso", subfmt="date_hm"),
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
                data.append(Table.read(surveyfile, path="data"))
                log.info(
                    "Read survey '{}' from {}: found {} pulsars".format(
                        survey, surveyfile, len(data[-1])
                    )
                )
            else:
                log.info(
                    "Cannot find cached survey '{}' at {}: loading...".format(
                        survey, surveyfile
                    )
                )
                try:
                    out = PulsarSurvey.read(
                        survey_name=survey,
                        survey_specs=Surveys[survey],
                    )
                except IndexError:
                    # when some pages are replaced this can happen
                    continue
                if out.data is None:
                    continue

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
        if not isinstance(self.data["P"], MaskedColumn):
            self.data["P"] = MaskedColumn(self.data["P"], mask=np.isnan(self.data["P"]))
        if not isinstance(self.data["DM"], MaskedColumn):
            self.data["DM"] = MaskedColumn(
                self.data["DM"], mask=np.isnan(self.data["DM"])
            )

    def search(
        self,
        coord: SkyCoord,
        radius: u.quantity.Quantity = 1 * u.deg,
        DM: float = None,
        DM_tolerance: float = 10,
        return_json: bool = False,
        return_native: bool = True,
        deduplicate: typing.Union[bool, str] = False,
    ):
        """
        returns an astropy.table object with the sources that match the given criteria
        or a JSON object if return_json=True

        Args:
            coord (SkyCoord): position to search around (icrs or galactic)
            radius (u.quantity.Quantity, optional): radius to search
            DM (float, optional): central DM constraint (if None, do no DM cut)
            DM_tolerance (float, optional): DM tolerance for constraint
            return_json (bool, optional): return JSON instead of Table
            return_native (bool, optional): return in the same coordinates as the input
            deduplicate (bool or str, optional): deduplicate results (if "hide", will suppress duplicates)

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
        dedup_select = np.arange(len(output))
        if deduplicate:
            deduplicate_table(output)
            if isinstance(deduplicate, str) and deduplicate.lower() == "hide":
                orig_length = len(output)
                dedup_select = output["Duplicate?"] == None
                output = output[dedup_select]

                output.remove_column("Duplicate?")
                log.debug(
                    "Deduplication removed {} pulsars".format(orig_length - len(output))
                )

        if coord.name == "galactic" and return_native:
            g = self.coord[i][good][dedup_select].galactic
            output.add_column(Column(g.l, name="l", format=".6f"), index=1)
            output.add_column(Column(g.b, name="b", format=".6f"), index=2)
            output.remove_columns(["RA", "Dec"])
        if not return_json:
            return output
        if np.any(output["DM"].mask):
            DM_output = output["DM"].data
            DM_output[output["DM"].mask] = -999
            output["DM"] = DM_output
        if np.any(output["P"].mask):
            P_output = output["P"].data
            P_output[output["P"].mask] = -999
            output["P"] = P_output
        output_dict = {}
        if coord.name == "galactic" and return_native:
            output_dict["searchl"] = {
                "display_name": "Search l (deg)",
                "value": coord.l.deg,
            }
            output_dict["searchb"] = {
                "display_name": "Search b (deg)",
                "value": coord.b.deg,
            }
            output_dict["searchcoord"] = {
                "display_name": "Search Coord",
                "value": coord.to_string(),
            }
        else:
            # return ICRS by default
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
            if coord.name == "galactic" and return_native:
                lon_key = "l"
                lon_value = {"display_name": "l (deg)", "value": row["l"]}
                lat_key = "b"
                lat_value = {"display_name": "b (deg)", "value": row["b"]}

            else:
                lon_key = "ra"
                lon_value = {"display_name": "RA (deg)", "value": row["RA"]}
                lat_key = "dec"
                lat_value = {"display_name": "Dec (deg)", "value": row["Dec"]}

            output_dict[key] = {
                lon_key: lon_value,
                lat_key: lat_value,
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
                    "value": row["retrieval date"].to_value("iso", subfmt="date_hm"),
                },
            }
            if deduplicate:
                if isinstance(deduplicate, bool) or (
                    isinstance(deduplicate, str) and deduplicate.lower() != "hide"
                ):
                    output_dict[key]["duplicate"] = {
                        "display_name": "Duplicate?",
                        "value": row["Duplicate?"],
                    }

        return output_dict
