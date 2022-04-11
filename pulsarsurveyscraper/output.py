import os
import numpy as np
import pdf2image
from astropy.table import Table
from astropy.time import Time
from fpdf import FPDF  # Library to create pdfs

import pulsarsurveyscraper


def make_output(
    result,
    format,
    query_text,
    query_url,
    directory=None,
    prefix="scraper",
    font="Helvetica",
):
    # Add pdf page settings and font styles
    pdf = FPDF("L", "mm", "A4")
    pdf.add_page()
    pdf.set_font(font, "B")
    pdf.set_auto_page_break(auto=True, margin=20)

    # Define cell widths
    cellw = (pdf.w - pdf.l_margin) / 10
    height = pdf.font_size

    # Set the title
    now = Time.now()
    pdf.cell(
        0,
        height,
        f"Scraper results for query submitted at {now.isot} (UTC)",
        align="C",
    )
    pdf.ln(height)
    pdf.set_font(font)  # Reset the font from bold to normal

    # Add the search query to the pdf summary
    pdf.multi_cell(0, 15, query_text, ln=1)

    # Append the result to the pdf file
    pdf.cell(0, height, "Found {} matches:".format(len(result)), ln=1)
    pdf.ln(height)

    # Write the data to the pdf file
    # But before that make some changes to the table data format
    # Create a dummy table to copy over the result and make changes
    pdf_res = Table.copy(result, copy_data=True)
    if "RA" in result.keys():
        pdf_res["RA"] = np.round(pdf_res["RA"], 6)
        pdf_res["Dec"] = np.round(pdf_res["Dec"], 6)
    else:
        pdf_res["l"] = np.round(pdf_res["l"], 6)
        pdf_res["b"] = np.round(pdf_res["b"], 6)
        
    pdf_res["P"] = np.round(pdf_res["P"], 6)
    pdf_res["Distance"] = np.round(pdf_res["Distance"], 6)

    units = np.array(["", "deg", "deg", "ms", "pc / cm^3", "", "", "deg"])

    # Set the widths of the columns
    width = np.ones(len(pdf_res.colnames)) * cellw
    width[5] *= 1.4
    width[6] *= 2.25

    # Write the header to the table
    pdf.set_font(font, "B")
    for c in range(2):
        # Write the column names first and then their units
        header = np.array([pdf_res.colnames, units])
        for i, head in enumerate(header[c]):
            pdf.cell(width[i], 2 * height, str(head), border=1, align="C")
        pdf.ln(2 * height)

    # Write the table data
    pdf.set_font(font)
    for row in pdf_res:
        for i, entry in enumerate(row):
            if i == 5:
                # This is the survey name entry solt where the link is to be added
                url = pulsarsurveyscraper.Surveys[entry]["url"]
                pdf.set_text_color(0, 0, 255)
                pdf.cell(
                    width[i], 2 * height, str(entry), border=1, link=url, align="C"
                )
                pdf.set_text_color(0, 0, 0)
            else:
                pdf.cell(width[i], 2 * height, str(entry), border=1, align="C")
        pdf.ln(2 * height)
    pdf.ln(2 * height)

    # Add it to the summary pdf
    pdf.cell(0, height, "API link:", ln=1)
    pdf.ln(height)
    pdf.set_text_color(0, 0, 255)
    pdf.cell(0, height, str(query_url), link=query_url, ln=1)

    if directory is None:
        directory = os.curdir

    timestamp = (f"{now.isot}").replace(":","_")
    if format == "pdf":
        outfile = os.path.join(directory, f"{prefix}_{timestamp}.pdf")
        pdf.output(outfile)
        return outfile
    elif format == "png":
        outfile = os.path.join(directory, f"{prefix}_{timestamp}.png")
        bytes = pdf.output()
        img = pdf2image.convert_from_bytes(bytes)
        if len(img) > 1:
            print(
                "PNG output only contains first page of results...please consider storing output in PDF format"
            )
        img[0].save(outfile)

        return outfile
