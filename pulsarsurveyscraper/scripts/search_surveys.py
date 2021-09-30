#!/usr/bin/env python
import sys
import os
import argparse
import json
import requests
import re
from astropy.coordinates import SkyCoord
from astropy import units as u
import pulsarsurveyscraper
from fpdf import FPDF   # Library to create pdfs
import numpy as np
from astropy.time import Time
from astropy.table import Table
import pdf2image


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
        "-p", "--pdf",
        action='store_true',
        default=True,
        help="Save the results as pdf file?"
    )
    parser.add_argument(
        "-pn", "--png",
        action='store_true',
        default=False,
        help="Save the results as png file?"
    )

    args = parser.parse_args()
    if args.verbosity == 1:
        pulsarsurveyscraper.log.setLevel("INFO")
    elif args.verbosity >= 2:
        pulsarsurveyscraper.log.setLevel("DEBUG")

    # Snippet for creating and rendering the pdf file

    # Add pdf page settings and font styles 
    pdf = FPDF('L', 'mm', 'A4')
    pdf.add_page()
    pdf.set_font('Helvetica', "B")
    pdf.set_auto_page_break(auto=True, margin=20)

    # Define cell widths
    cellw = (pdf.w - pdf.l_margin)/10
    height = pdf.font_size

    # Set the title
    now = Time.now()
    pdf.cell(0, height, f'Summary of survey results for the query submitted at {now.isot} (UTC)', align="C")
    pdf.ln(height)
    pdf.set_font('Helvetica') # Reset the font from bold to normal

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
    if not args.galactic:
        search_query_txt = "Searching {:.1f}deg around RA,Dec = {} = {}d,{}d".format(
                args.radius,
                coord.to_string("hmsdms", sep=":"),
                coord.ra.to_string(decimal=True),
                coord.dec.to_string(decimal=True, alwayssign=True),
            )

    else:
        search_query_txt = "Searching {:.1f}deg around l,b = {}d,{}d".format(
                args.radius,
                coord.l.to_string(decimal=True),
                coord.b.to_string(decimal=True, alwayssign=True),
            )

    print(search_query_txt)
    # Add the same to the pdf summary
    pdf.multi_cell(0, 15, search_query_txt, ln=1)

    if args.dm is not None and args.dmtol is not None:
        print(
            "Also requiring DM with +/-{:.1f} of {:.1f} pc/cm**2".format(
                args.dmtol, args.dm
            )
        )
    result = pulsar_table.search(
        coord,
        radius=args.radius * u.deg,
        DM=args.dm,
        DM_tolerance=args.dmtol,
        return_json=args.json,
        deduplicate=args.dedup,
    )

    # Append the result to the pdf file
    pdf.cell(0, height, "Found {} matches:".format(len(result)), ln=1)
    pdf.ln(height)

    if not args.json:
        print("Found {} matches:".format(len(result)))
    else:
        print("Found {} matches:".format(result["nmatches"]))
        result = json.dumps(result)
    print(result)

    # Write the data to the pdf file
    # But before that make some changes to the table data format 
    # Create a dummy table to copy over the result and make changes
    pdf_res = Table.copy(result, copy_data=True)
    pdf_res['RA'] = np.round(pdf_res['RA'], 6)
    pdf_res['Dec'] = np.round(pdf_res['Dec'], 6)
    pdf_res['P'] = np.round(pdf_res['P'], 6)
    pdf_res['Distance'] = np.round(pdf_res['Distance'], 6)

    units = np.array(['', 'deg', 'deg', 'ms', 'pc / cm^3', '', '', 'kpc'])

    # Set the widths of the columns
    width = np.ones(len(pdf_res.colnames))*cellw
    width[5] *=1.4
    width[6] *=2.25

    # Write the header to the table
    pdf.set_font('Helvetica', 'B')
    for c in range(2):
        # Write the column names first and then their units
        header = np.array([pdf_res.colnames, units])
        for i, head in enumerate(header[c]):
            pdf.cell(width[i], 2*height, str(head), border=1, align='C')
        pdf.ln(2*height)    


    # Write the table data
    pdf.set_font('Helvetica')
    for row in pdf_res:
        for i, entry in enumerate(row):
            if i==5: 
                # This is the survey name entry solt where the link is to be added
                url = pulsarsurveyscraper.Surveys[entry]['url']
                pdf.set_text_color(0,0,255)
                pdf.cell(width[i], 2*height, str(entry), border=1, link=url, align='C')
                pdf.set_text_color(0,0,0)
            else:
                pdf.cell(width[i], 2*height, str(entry), border=1, align='C')
        pdf.ln(2*height)
    pdf.ln(2*height)

    # Create the api link to be appended at the end of PDF for 
    # reproducibility
    if args.dm is not None:
        query_dict = {'type':'search', 'ra':str(ra), 'dec':str(dec), 
            'radius':args.radius, 'dm':args.dm, 'dmtol':args.dmtol}
    elif args.dm is None:
        query_dict = {'type':'search', 'ra':str(ra), 'dec':str(dec), 
            'radius':args.radius}

    # Get the global URL for the query web page:
    query_url = "https://pulsar.cgca-hub.org/api"
    response = requests.get(query_url, params=query_dict)

    # Add it to the summary pdf
    pdf.cell(0, height, "API link:", ln=1)
    pdf.ln(height)
    pdf.set_text_color(0, 0, 255)
    pdf.cell(0, height, str(response.url), link=response.url, ln=1)
    
    if not args.png:
        outfile=f'summary_{now.isot}.pdf'
        pdf.output(outfile)
    elif args.png:
        bytes = pdf.output(dest='S')
        img = pdf2image.convert_from_bytes(bytes)
        if len(img)>1:
            print("Search query returned multiple sources that have to be stored in multiple pages/files and hence just the first page/file is being stored as image file... Please consider storing it in pdf format")
        img[0].save(f'summary_{now.isot}.png')


if __name__ == "__main__":
    main()
