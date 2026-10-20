from data_analysis import compute_all, format_result
from exact_diagonalization import block_H, find_all_observables
from make_jobs import make_file_name
from pathlib import Path

# By ChatGPT:
###
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor
import numpy as np

def add_results_table(doc, table_number, beta, J_1, J_2,
                      ed, means_z, errs_z, means_x, errs_x):

    # Table caption
    p = doc.add_paragraph()
    p.add_run(
        f"Table {table_number}. "
        f"β = {beta}, J₁ = {J_1}, J₂ = {J_2}."
    )

    # Your sample has 15 observables + header
    rows = [
        ("E",    "E",     "E",     "E"),
        ("M_z²",  "mz2",   "M^2",   None),
        ("M_z⁴",  "mz4",   "M^4",   None),
        ("M_x²",  "mx2",   None,    "M^2"),
        ("M_x⁴",  "mx4",   None,    "M^4"),
        ("C",    "C",     "C",     "C"),
        ("χz",   "chi_z", "chi_M", None),
        ("χx",   "chi_x", None,    "chi_M"),
        ("χj",   "chi_j", "chi_J", "chi_J"),
        ("ξz",   "xi_z",  "xi_M",  None),
        ("ξx",   "xi_x",  None,    "xi_M"),
        ("A",    "A",     "A",     "A"),
        ("B",    "B",     "B",     "B"),
        ("ZZ",   "ZZ",    "ZZ",    "ZZ"),
        ("XX",   "XX",    "XX",    "XX"),
    ]

    table = doc.add_table(rows=1, cols=4)
    table.style = "Table Grid"
    table.autofit = True

    # Header
    headers = ["Observable", "ED", "SSE (Z)", "SSE (X)"]

    for j, header in enumerate(headers):
        set_cell(table.rows[0].cells[j], header, bold=True)

    # Rows
    for label, ed_key, z_key, x_key in rows:
        cells = table.add_row().cells

        # Observable
        set_cell(cells[0], label, align=WD_ALIGN_PARAGRAPH.LEFT)

        # ED
        if ed_key is None:
            set_cell(cells[1], "-")
        else:
            set_cell(cells[1], str(ed[ed_key]))

        # SSE Z
        if z_key is None:
            set_cell(cells[2], "-")
        else:
            mean_z = means_z[z_key]
            err_z = errs_z[z_key]
            z_text = format_result(mean_z, err_z)

            # Red if ED is outside the SSE ± error interval
            z_color = None
            if ed_key is not None and np.isfinite(mean_z) and np.isfinite(err_z):
                if not (mean_z - 3 * err_z <= ed[ed_key] <= mean_z + 3 * err_z):
                    z_color = RGBColor(255, 0, 0)

            set_cell(
                cells[2],
                z_text,
                color=z_color
            )

        # SSE X
        if x_key is None:
            set_cell(cells[3], "-")
        else:
            mean_x = means_x[x_key]
            err_x = errs_x[x_key]
            x_text = format_result(mean_x, err_x)

            # Red if ED is outside the SSE ± error interval
            x_color = None
            if ed_key is not None and np.isfinite(mean_x) and np.isfinite(err_x):
                if not (mean_x - 3 * err_x <= ed[ed_key] <= mean_x + 3 * err_x):
                    x_color = RGBColor(255, 0, 0)

            set_cell(
                cells[3],
                x_text,
                color=x_color
            )

    # Blank paragraph after the table
    doc.add_paragraph()

def set_cell(cell, text, bold=False, align=WD_ALIGN_PARAGRAPH.CENTER, color=None):
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.alignment = align

    run = paragraph.add_run(str(text))
    run.bold = bold
    run.font.size = Pt(10)

    if color is not None:
        run.font.color.rgb = color

"""
doc = Document()          # <-- Create a new file
table_number = 1
###

# Pre-compute ED
L = 2
op_blks, mz_blks, mx_blks, sectors = block_H(L)

# Loop through all parameter points
for beta in [0.1, 1, 5]:
    for J_1 in [0.1, 1, 5]:
        for J_2 in [0, 0.1, 1, 5]:

            if J_2 > J_1:
                continue
            
            print(f"Starting: L=2, beta={beta} J1={J_1}, J2={J_2}")
            # Compute ED results
            ed_results = find_all_observables(op_blks, mz_blks, mx_blks, sectors, L, J_1, J_2, beta)

            # Compute SSE results from Z-basis
            filename_Z = make_file_name(L, beta, J_1, J_2, "Z", 0, "simulate") + "0.h5"
            means_z, errs_z = compute_all(filename_Z)

            # Compute SSE results from X-basis
            filename_X = make_file_name(L, beta, J_1, J_2, "X", 0, "simulate") + "0.h5"
            means_x, errs_x = compute_all(filename_X)

            # Make docx table
            # By ChatGPT
            ###
            add_results_table(
                doc,
                table_number,
                beta,
                J_1,
                J_2,
                ed_results,
                means_z,
                errs_z,
                means_x,
                errs_x,
            )

            table_number += 1
            ###

doc.save(Path("Analysis") / "SSE_ED_results.docx")
"""

# Compute SSE results from Z-basis
filename_Z = make_file_name(3, 1, 1, 1, "Z", 0, "simulate") + "1.h5"
means_z, errs_z = compute_all(filename_Z)
for observable in means_z.keys():
    print(observable, format_result(means_z[observable], errs_z[observable]))