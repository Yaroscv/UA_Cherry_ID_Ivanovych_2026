#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03.UA_Cherry_ID_Harmonise_ECPGR  —  harmonise Ukrainian SSR calls to the ECPGR allele scale
==========================================================================================

Project : UA sweet cherry ID (Y. Ivanovych) — manuscript §2.3, Supplementary Table S3
Version : see VERSION below (semantic versioning; bump on every edit, log it in CHANGELOG)
Author  : Yaroslav Ivanovych; code consolidated with AI assistance (Claude) from the
          as-run steps of 2026-07-18 archived in etc/03.UA_Cherry_ID_Harmonise_asrun_2026-07-18/

What it does
------------
1. Reads the raw Ukrainian genotypes (94 accessions, 18 SSR + 3 S-locus; only the six loci
   shared with the ECPGR aligned set are used: BPPCT034, BPPCT037, CPSCT038, EMPa002,
   EMPaS02, EMPaS14).
2. Builds a per-allele correspondence UKR raw bp -> ECPGR bp from the anchor cultivars typed
   both here and in the European sets (rank pairing per locus, consensus = majority vote).
3. Maps every Ukrainian allele: anchor value if the allele is an anchor; linear interpolation
   between anchors; extrapolation (offset of the nearest end anchor) outside the span; then
   snaps to the nearest allele of the COMBINED ECPGR ladder (ties -> the lower allele).
4. Repeats step 3 against the narrower Barreneche-314 ladder and reports the differences.
5. Flags alleles: "Extrapolated" (outside the calibration span) or "Ambiguous" (interpolated
   value midway, within 0.5 bp, between two ladder alleles). Interpolations are not flagged.
6. Writes Table S3 (xlsx, same content as Supplementary_Table_..._harmonisation_v2.xlsx of
   2026-07-18), the 94 x 6 harmonised matrix (txt), diagnostics (tsv) and a run manifest
   (json: version, input/output SHA-256, counts, check results).
7. Regression checks (default on): the run must reproduce the counts of Table S3 v2.

Usage
-----
    python 03.UA_Cherry_ID_Harmonise_ECPGR.py                 # ROOT = folder of this script
    python 03.UA_Cherry_ID_Harmonise_ECPGR.py --root C:/R/UA_Cherry_ID --out results/harmonise
    python 03.UA_Cherry_ID_Harmonise_ECPGR.py --no-check
    python tests/test_03_harmonise.py                          # regression test vs Table S3 v2

Versioning: the version lives in VERSION (printed, and written to manifest.json), not in the
file name. MAJOR = method change (results may differ); MINOR = new options/outputs, results
unchanged; PATCH = fixes and documentation. Superseded versions: etc/03.UA_Cherry_ID_Harmonise_ECPGR_v<X.Y.Z>.py

Requires Python >= 3.8 and openpyxl.

CHANGELOG
---------
1.0.1  2026-11-12  File renamed by the author to 03.UA_Cherry_ID_Harmonise_ECPGR.py (version kept
                   in VERSION, not in the name); usage text updated; versioning rule stated.
                   No code change; the test still reproduces Table S3 v2 cell for cell.
1.0.0  2026-11-12  First versioned release. Consolidates the six as-run steps of 2026-07-18
                   (recal_full.py, the recalibrated-txt writer, the flag rule, the Table S3
                   builder and its two patches) into one script; no hard-coded results.
                   Changes against the as-run code, each tested to leave every output cell
                   unchanged:
                   (a) reference names are stripped on BOTH sides of the lookup — as run,
                       'FRA Regina ' (trailing space in the query) never matched, so Regina
                       was silently not used as an anchor (13 accessions used, not 14);
                   (b) paths are relative to --root; inputs are hashed into the manifest;
                   (c) the ladder-diff and anchor-verification tables are written as tsv.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import platform
import sys
from collections import Counter, defaultdict
from pathlib import Path

VERSION = "1.0.1"
SCRIPT = "03.UA_Cherry_ID_Harmonise_ECPGR"

# ── Inputs (relative to <root>/data) ────────────────────────────────────────────────────
F_UKR   = "UKR_Sw_Cherry_ID_18_SSR_3_S-locus_data.txt"          # raw UKR matrix (2022-05-10)
F_B314  = "EU Cherry Barreneche, 2021, 14 SSR - 314 genotypes. 14 SSR data.txt"
F_ORD47 = "ECPGR. Sw Cherry Ordidge, 2021. 14 SSR data.txt"     # author's 47-accession subset of Ordidge 2021
F_LCOMB = "AFP_ECPGR_Sw.Cherry ID - Ordidge_Barreneche, 2021, ECPGR - Familias_recalc.txt"  # combined ladder (2,394)
F_LB314 = "AFP_EU Cherry Barreneche, 2021, SSR. 14 SSR allele freq. list.txt"               # Barreneche-314 ladder

# ── Outputs (in <out>) ──────────────────────────────────────────────────────────────────
O_XLSX  = "Supplementary_Table_S3_UKR_SSR_ECPGR_harmonisation.xlsx"
O_TXT   = "UKR_Cherry_94_harmonised_ECPGR.txt"
O_DIFF  = "ladder_diff_Barreneche314_vs_combined.tsv"
O_ANCH  = "anchor_verification.tsv"
O_MAP   = "allele_correspondence.tsv"
O_MANI  = "manifest.json"

# Loci: spelling in the raw UKR file -> spelling in the ECPGR files
LOCI = [("BPPCT034", "BPPCT034"), ("BPPCT037", "BPPCT037"), ("CPSCT038", "CPSCT038"),
        ("EMPA002", "EMPa002"), ("EMPaS02", "EMPaS02"), ("EMPaS14", "EMPaS14")]

# Anchor cultivars: (UKR sample name, Barreneche-314 name, Ordidge-47 name).
# Barreneche is used when named; otherwise Ordidge-47. ORDER MATTERS: it fixes the
# insertion order of the votes and hence the tie-break of the consensus.
ANCHORS = [
    ("Napoleon 92", "FRA Napoleon", "Napoleon-STD"),
    ("Xapata 130", "FRA Xapata", None),
    ("Jaboulay 108", None, "Jaboulay-France"),
    ("Regina REG", "FRA Regina", "Regina-France"),
    ("Burlat 77", "FRA Bigarreau Hatif Burlat", "Bigarreau Hatif Burlat-France"),
    ("Garnet 97", "FRA Garnet", "Garnet-France"),
    ("Lapins LAP", None, "Lapins-Germany"),
    ("Saint Georges 107", "FRA Saint Georges", "Saint Georges-France"),
    ("Cristobalina 19", "FRA Cristobalina", "Cristobalina-France"),
    ("Folfer® 223", None, "Folfer-France"),
    ("Adriana 220", "ITA Adriana", None),
    ("Dönissens Gelbe 110", "DEU Drogans Gelbe Knorpelkirsche", "Drogans Gelbe Knorpel-Germany"),
    ("Drogan's Yellow DRO/1", None, "Donissens Gelbe Knorpelkirsche_SG 45-Germany"),
    ("Margit 10", None, "Margit-Germany"),
]

# Regression targets = Table S3 v2 (2026-07-18). A run that differs fails the check.
EXPECTED = {
    "n_accessions": 94,
    "n_anchors_used": 14,                      # 13 as run (Regina bug); map identical
    "map_sizes": {"BPPCT034": 8, "BPPCT037": 7, "CPSCT038": 5, "EMPa002": 2, "EMPaS02": 7, "EMPaS14": 3},
    "spans": {"BPPCT034": [232, 274], "BPPCT037": [154, 180], "CPSCT038": [202, 222],
              "EMPa002": [124, 126], "EMPaS02": [152, 166], "EMPaS14": [218, 232]},
    "ladder_sizes_combined": {"BPPCT034": 24, "BPPCT037": 25, "CPSCT038": 7, "EMPa002": 19, "EMPaS02": 25, "EMPaS14": 19},
    "ladder_sizes_b314": {"BPPCT034": 15, "BPPCT037": 10, "CPSCT038": 6, "EMPa002": 5, "EMPaS02": 8, "EMPaS14": 7},
    "ladder_diff_accessions": 4,
    "ladder_diff_cells": 8,
    "flagged_accessions": 8,
    "flagged_alleles": 18,
}

FLAG_EXTRAP = "Extrapolated (outside calibration range)"
FLAG_AMBIG = "Ambiguous (midway between two ECPGR alleles)"


# ════════════════════════════════════════════════════════════════════════════════════════
# Reading
# ════════════════════════════════════════════════════════════════════════════════════════
def _int(x) -> int:
    try:
        return int(x)
    except (TypeError, ValueError):
        return 0


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def read_ukr(path: Path) -> tuple[list[str], dict]:
    """Raw UKR matrix -> (sample order, {sample: {raw_locus: (a1, a2)}}); '*' or blank -> 0."""
    lines = read_lines(path)
    hdr = lines[0].split("\t")
    col = {raw: hdr.index(raw + " 1") for raw, _ in LOCI}
    order, geno = [], {}
    for ln in lines[1:]:
        p = ln.split("\t")
        if not p or not p[0]:
            continue
        nm = p[0]
        order.append(nm)
        geno[nm] = {raw: (_int(p[col[raw]]), _int(p[col[raw] + 1])) for raw, _ in LOCI}
    return order, geno


def read_ref(path: Path) -> dict:
    """Reference genotype table -> {stripped name: {raw_locus: (a1, a2)}}."""
    rows = [l.split("\t") for l in read_lines(path)]
    hdr = rows[0]
    out = {}
    for r in rows[1:]:
        if not r or not r[0].strip():
            continue
        g = {}
        for raw, ec in LOCI:
            i = hdr.index(ec + " 1")
            g[raw] = (_int(r[i].strip()), _int(r[i + 1].strip()))
        out.setdefault(r[0].strip(), g)          # first occurrence wins (as run)
    return out


def read_ladder(path: Path) -> dict[str, list[int]]:
    """Familias frequency list -> {ECPGR locus: [alleles in file order]}."""
    d, cur = {}, None
    for line in read_lines(path):
        s = line.split("\t")
        if s[0] and (len(s) < 2 or s[1] == ""):
            cur = s[0].strip().replace("EMPA002", "EMPa002")
            d[cur] = []
            continue
        if cur and len(s) >= 2 and s[0] and s[1]:
            try:
                d[cur].append(int(float(s[0])))
            except ValueError:
                pass
    return d


# ════════════════════════════════════════════════════════════════════════════════════════
# Calibration
# ════════════════════════════════════════════════════════════════════════════════════════
def anchor_reference(name_b, name_o, b314, ord47):
    if name_b:                                   # as run: Barreneche if named, else Ordidge-47
        return b314.get(name_b.strip())
    if name_o:
        return ord47.get(name_o.strip())
    return None


def build_map(ukr, b314, ord47):
    votes = {raw: defaultdict(Counter) for raw, _ in LOCI}
    used = []
    for un, bn, on in ANCHORS:
        if un not in ukr:
            continue
        ref = anchor_reference(bn, on, b314, ord47)
        if not ref:
            continue
        used.append((un, bn or on, "Barreneche-314" if bn else "Ordidge-47"))
        for raw, _ in LOCI:
            u = sorted(a for a in ukr[un][raw] if a)
            r = sorted(a for a in ref[raw] if a)
            if len(u) == len(r) and u:
                for ua, ra in zip(u, r):
                    votes[raw][ua][ra] += 1
    amap = {raw: {ua: cnt.most_common(1)[0][0] for ua, cnt in votes[raw].items()} for raw, _ in LOCI}
    return amap, used


def interp(amap_l: dict, a: int):
    """-> (EU-scale value, kind) with kind in {'anchor', 'interp', 'extrap'}."""
    if a in amap_l:
        return amap_l[a], "anchor"
    an = sorted(amap_l.items())
    us = [u for u, _ in an]
    es = [e for _, e in an]
    if not us:
        return a, "anchor"
    if a < us[0]:
        return es[0] + (a - us[0]), "extrap"
    if a > us[-1]:
        return es[-1] + (a - us[-1]), "extrap"
    for i in range(len(us) - 1):
        if us[i] < a < us[i + 1]:
            return es[i] + (a - us[i]) * (es[i + 1] - es[i]) / (us[i + 1] - us[i]), "interp"
    return a, "anchor"


def snap(ladder_sorted: list[int], eu: float) -> int:
    return min(ladder_sorted, key=lambda x: abs(x - eu))     # ties -> lower allele


def harmonise(ukr, amap, ladders):
    """-> {ladder name: {sample: {raw_locus: (h1, h2)}}}; 0 stays 0 (missing)."""
    out = {}
    for lname, lad in ladders.items():
        out[lname] = {}
        for nm, g in ukr.items():
            out[lname][nm] = {}
            for raw, ec in LOCI:
                ls = sorted(set(lad[ec]))
                out[lname][nm][raw] = tuple(snap(ls, interp(amap[raw], a)[0]) if a else 0 for a in g[raw])
    return out


def flag_alleles(order, ukr, amap, ladder_comb):
    rows, flagged = [], []
    for nm in order:
        for raw, ec in LOCI:
            lad = sorted(ladder_comb[ec])
            for a in ukr[nm][raw]:
                if not a:
                    continue
                eu, kind = interp(amap[raw], a)
                ds = sorted(lad, key=lambda x: abs(x - eu))
                snapped = ds[0]
                tie = len(ds) > 1 and abs(abs(ds[0] - eu) - abs(ds[1] - eu)) < 0.5
                flag = FLAG_EXTRAP if kind == "extrap" else (FLAG_AMBIG if kind == "interp" and tie else None)
                if flag:
                    rows.append([nm, ec, a, snapped, flag])
                    if nm not in flagged:
                        flagged.append(nm)
    return rows, flagged


# ════════════════════════════════════════════════════════════════════════════════════════
# Table S3 workbook (content and layout of Table S3 v2, 2026-07-18)
# ════════════════════════════════════════════════════════════════════════════════════════
NOTES = [  # (text, size, bold) — frozen text of Table S3 v2; edit = new VERSION
    ("Supplementary Table S_. Harmonisation of Ukrainian sweet cherry (Prunus avium L.) SSR genotypes", 12, True),
    ("to the combined ECPGR reference allele scale.", 12, True), ("", 10, False),
    ("Source: raw Ukrainian genotypes (18 SSR); the three S-locus markers (PaConsI, PaConsII, F-box) were", 10, False),
    ("removed. Six loci overlap the ECPGR standard cherry set (Clarke & Tobutt, 2009) and are reported", 10, False),
    ("here: BPPCT034, BPPCT037, CPSCT038, EMPa002, EMPaS02, EMPaS14.", 10, False), ("", 10, False),
    ("Reference ladder: the COMBINED ECPGR set (Ordidge et al., 2021 + Barreneche et al., 2021; 2,394", 10, False),
    ("genotypes), wider than the Barreneche-only 314-genotype set (24/25/7/19/25/19 vs 15/10/6/5/8/7", 10, False),
    ("alleles at these six loci). A per-allele correspondence (sheet \"Allele_correspondence\") was built", 10, False),
    ("from 13 cultivars genotyped both here and in the European collections; Ukrainian alleles were mapped", 10, False),
    ("to the nearest allele of the combined ladder.", 10, False), ("", 10, False),
    ("Recalibration was run for all 94 accessions directly from the raw data. Only 4 accessions differ", 10, True),
    ("from a Barreneche-314-based calibration (Heneralska, Plena, AVCV4, AVCV5; 8 locus-genotypes); the", 10, False),
    ("remaining 90 are identical on both ladders. Differences arise only for alleles lying outside the", 10, False),
    ("range spanned by the calibration cultivars.", 10, False), ("", 10, False),
    ("All accessions are RETAINED. Sheet \"Flagged_alleles\" lists individual alleles whose harmonised value", 10, False),
    ("is less certain. IMPORTANT: \"Extrapolated\" means the Ukrainian allele lies OUTSIDE the span of", 10, True),
    ("Ukrainian alleles seen in the shared calibration cultivars (the \"calibration span\", given per locus", 10, False),
    ("in sheet \"Allele_correspondence\"), so the correspondence had to be extended beyond the calibrated", 10, False),
    ("range. It does NOT mean the resulting allele is absent from the ECPGR ladder - the harmonised value", 10, False),
    ("is usually a perfectly valid ECPGR allele (e.g. EMPaS02 raw 168 -> 151: the span ends at UKR 166,", 10, False),
    ("yet 151 exists on the ladder). The flag reflects calibration coverage, not ladder membership.", 10, False),
    ("\"Ambiguous\" means the mapped value falls midway between two ECPGR alleles. Alleles interpolated", 10, False),
    ("BETWEEN calibration anchors are NOT flagged. This is a data-quality annotation, not an exclusion.", 10, False),
    (None, 11, False),
    ("Note on the pooled reference ladder", 10, True),
    ("The combined ECPGR ladder pools sub-datasets aligned by Ordidge et al. (2021). Those authors report that", 10, False),
    ("residual differences of 1-2 bp remain after alignment and \"all of those that were accepted were", 10, False),
    ("conceivably caused by either rounding or calling error\"; that for seven loci the alignment factor was", 10, False),
    ("\"dependent on whether alleles were scored as odd or even\"; that EMPaS02 in particular required an", 10, False),
    ("allele-range-specific factor (alleles of 146 and 148 bp adjusted by 4 bp, the remainder by 3 bp); and that", 10, False),
    ("24 queries remained unresolved. They further caution that \"any increase in allele number that cannot be", 10, False),
    ("accounted for might be investigated further to confirm whether it could be an artefact of either", 10, False),
    ("inconsistent scoring or alignment\", and that data aligned to the published reference sizes \"will require a", 10, False),
    ("further stage of realignment\" to their dataset.", 10, False), ("", 10, False),
    ("Consistent with this, the pooled ladder contains adjacent 1 bp entries at several loci (EMPaS02 has the", 10, False),
    ("most: 15 adjacent pairs). At EMPaS02, 151, 152 and 153 all appear, but 152 occurs in only three records", 10, False),
    ("(all from the Barreneche et al. 2021 profiles merged here) and 153 in a single record (Italian national", 10, False),
    ("data); no accession carries both, so they cannot be tested against each other and may represent the same", 10, False),
    ("allele scored on two sub-datasets that were never cross-aligned. The pooled ladder should therefore not be", 10, False),
    ("read as a clean biological allele series. This does not affect the present harmonisation: no Ukrainian", 10, False),
    ("allele maps above 151 at EMPaS02 (harmonised values span 134-151 only).", 10, False),
]

FLAG_SHEET_NOTE = ('Individual alleles with less certain harmonisation. All accessions are RETAINED. '
                   '"Extrapolated" = the raw UKR allele lies outside the calibration span (last column), '
                   'NOT that the harmonised allele is missing from the ECPGR ladder.')
AC_FOOTNOTE = ('Note: the "Full combined ECPGR ladder" is a pooled list from sub-datasets aligned by Ordidge et al. (2021); '
               'adjacent 1 bp entries may reflect rounding or odd/even scoring rather than distinct alleles (see Notes).')
LADDER_CAVEAT = [
    [151, 'many', 'both sub-datasets', 'yes', 'yes', 'well supported; the extrapolated Ukrainian allele (raw 168) maps here'],
    [152, 3, 'Barreneche et al. 2021 profiles merged here (Zalanka, Bigarreau Camus de Venasque, Aelita)', 'yes', 'no', 'rare; sub-dataset specific'],
    [153, 1, 'Italian national data (Ravenna tardiva)', 'no', 'yes', 'rare; sub-dataset specific'],
    ['-', '-', '-', '-', '-', '152 and 153 never co-occur in one accession, so they cannot be tested against each other. '
     'They may be the same allele scored on two sub-datasets that were not cross-aligned (cf. Ordidge et al. 2021, who note '
     'residual rounding/odd-even effects and that a further stage of realignment is required).'],
]


def write_table_s3(path: Path, order, harm_comb, amap, ladder_comb, flag_rows, flagged, spans):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter

    ari = lambda **k: Font(name="Arial", **k)
    HDR = PatternFill("solid", fgColor="1F4E78"); HF = ari(bold=True, color="FFFFFF")
    SUB = PatternFill("solid", fgColor="DDEBF7"); YEL = PatternFill("solid", fgColor="FFF2CC")
    th = Side(style="thin", color="BFBFBF"); BORD = Border(th, th, th, th)

    def head(ws, row, cols):
        for j in cols:
            c = ws.cell(row=row, column=j)
            c.fill, c.font, c.border = HDR, HF, BORD
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

    wb = Workbook()
    # Notes
    s = wb.active; s.title = "Notes"
    for i, (t, z, b) in enumerate(NOTES, 1):
        s.cell(row=i, column=1, value=t if t else None).font = ari(size=z, bold=b)
    s.column_dimensions["A"].width = 104
    # Harmonised_genotypes
    s = wb.create_sheet("Harmonised_genotypes")
    top = ["No.", "Accession"]
    for _, ec in LOCI:
        top += [ec, ""]
    s.append(top + ["Note"])
    for i in range(6):
        c0 = 3 + 2 * i
        s.merge_cells(start_row=1, start_column=c0, end_row=1, end_column=c0 + 1)
    s.merge_cells(start_row=1, start_column=15, end_row=2, end_column=15)
    s.append(["", ""] + ["Allele 1", "Allele 2"] * 6 + [""])
    for n, nm in enumerate(order, 1):
        row = [n, nm]
        for raw, _ in LOCI:
            a = sorted(harm_comb[nm][raw])
            row += [a[0] if a[0] else "-", a[1] if a[1] else "-"]
        row.append("flagged allele(s) - see Flagged_alleles" if nm in flagged else "")
        s.append(row)
    head(s, 1, [1, 2, 3, 5, 7, 9, 11, 13, 15])
    for j in range(1, 15):
        c = s.cell(row=2, column=j); c.fill = SUB; c.font = ari(bold=True, size=9); c.border = BORD
        c.alignment = Alignment(horizontal="center")
    for r in range(3, s.max_row + 1):
        for j in range(1, 16):
            c = s.cell(row=r, column=j); c.font = ari(size=9); c.border = BORD
            if 2 < j < 15:
                c.alignment = Alignment(horizontal="center")
        if s.cell(row=r, column=15).value:
            for j in range(1, 16):
                s.cell(row=r, column=j).fill = YEL
    s.column_dimensions["A"].width = 5; s.column_dimensions["B"].width = 30
    for j in range(3, 15):
        s.column_dimensions[get_column_letter(j)].width = 8
    s.column_dimensions["O"].width = 38; s.freeze_panes = "C3"
    # Allele_correspondence
    s = wb.create_sheet("Allele_correspondence")
    s.append(["Locus", "Ukrainian allele (bp)", "ECPGR allele (bp)", "Full combined ECPGR ladder (bp)",
              "Calibration span (UKR bp)"])
    for raw, ec in LOCI:
        lad = ",".join(map(str, ladder_comb[ec]))
        for ua in sorted(amap[raw]):
            s.append([ec, ua, amap[raw][ua], lad, f"{spans[ec][0]}-{spans[ec][1]}"])
    head(s, 1, range(1, 6))
    for r in range(2, s.max_row + 1):
        for j in range(1, 6):
            c = s.cell(row=r, column=j); c.font = ari(size=9); c.border = BORD
            if j in (2, 3, 5):
                c.alignment = Alignment(horizontal="center")
    for col, w in zip("ABCDE", [12, 18, 16, 60, 20]):
        s.column_dimensions[col].width = w
    s.freeze_panes = "A2"
    s.cell(row=s.max_row + 2, column=1, value=AC_FOOTNOTE).font = ari(size=9, italic=True)
    # Flagged_alleles
    s = wb.create_sheet("Flagged_alleles")
    s.append([FLAG_SHEET_NOTE]); s.append([])
    s.append(["Accession", "Locus", "Raw allele (bp)", "Harmonised ECPGR allele (bp)", "Flag", "Calibration span (UKR bp)"])
    for r in flag_rows:
        s.append(r + [f"{spans[r[1]][0]}-{spans[r[1]][1]}"])
    s.cell(row=1, column=1).font = ari(bold=True, size=10, italic=True)
    head(s, 3, range(1, 7))
    for r in range(4, s.max_row + 1):
        for j in range(1, 7):
            c = s.cell(row=r, column=j); c.font = ari(size=9); c.border = BORD; c.fill = YEL
            if j in (3, 4, 6):
                c.alignment = Alignment(horizontal="center")
    for col, w in zip("ABCDEF", [24, 12, 16, 26, 44, 20]):
        s.column_dimensions[col].width = w
    s.freeze_panes = "A4"
    # Ladder_caveat
    s = wb.create_sheet("Ladder_caveat")
    s.append(["EMPaS02 upper-range alleles in the pooled ECPGR ladder - provenance and interpretation"]); s.append([])
    s.append(["Allele (bp)", "n records in combined DB", "Source of those records", "In Barreneche-314?",
              "In Ordidge national data?", "Interpretation"])
    for x in LADDER_CAVEAT:
        s.append(x)
    s.cell(row=1, column=1).font = ari(bold=True, size=11)
    head(s, 3, range(1, 7))
    for rz in range(4, s.max_row + 1):
        for j in range(1, 7):
            c = s.cell(row=rz, column=j); c.font = ari(size=9); c.border = BORD; c.fill = YEL
            c.alignment = Alignment(wrap_text=True, vertical="top", horizontal="center" if j in (1, 2) else "left")
    for col, w in zip("ABCDEF", [12, 20, 46, 18, 22, 62]):
        s.column_dimensions[col].width = w
    wb.save(path)


# ════════════════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════════════════
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(map(str, header)) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v) for v in r) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=f"{SCRIPT} v{VERSION}")
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent), help="project root (contains data/)")
    ap.add_argument("--data", default="data", help="input folder, relative to --root")
    ap.add_argument("--out", default="results/harmonise", help="output folder, relative to --root")
    ap.add_argument("--no-check", action="store_true", help="skip the regression checks")
    a = ap.parse_args(argv)
    root = Path(a.root); data = root / a.data; out = root / a.out
    out.mkdir(parents=True, exist_ok=True)
    print(f"{SCRIPT} v{VERSION}  |  python {platform.python_version()}  |  root {root}")

    inputs = {k: data / v for k, v in
              {"ukr": F_UKR, "b314": F_B314, "ord47": F_ORD47, "ladder_comb": F_LCOMB, "ladder_b314": F_LB314}.items()}
    for k, p in inputs.items():
        if not p.exists():
            sys.exit(f"missing input {k}: {p}")

    order, ukr = read_ukr(inputs["ukr"])
    b314, ord47 = read_ref(inputs["b314"]), read_ref(inputs["ord47"])
    lad_comb, lad_b314 = read_ladder(inputs["ladder_comb"]), read_ladder(inputs["ladder_b314"])

    amap, used = build_map(ukr, b314, ord47)
    spans = {ec: [min(amap[raw]), max(amap[raw])] for raw, ec in LOCI}
    H = harmonise(ukr, amap, {"COMB": lad_comb, "B314": lad_b314})
    flag_rows, flagged = flag_alleles(order, ukr, amap, lad_comb)

    diff = []
    for nm in order:
        for raw, ec in LOCI:
            b, c = tuple(sorted(H["B314"][nm][raw])), tuple(sorted(H["COMB"][nm][raw]))
            if b != c:
                diff.append([nm, ec, "/".join(map(str, ukr[nm][raw])), "/".join(map(str, b)), "/".join(map(str, c))])

    anchor_rows = []
    for un, refname, src in used:
        ref = b314.get(refname.strip()) if src == "Barreneche-314" else ord47.get(refname.strip())
        ok = tot = 0
        for raw, ec in LOCI:
            h = tuple(sorted(x for x in H["COMB"][un][raw] if x)); r = tuple(sorted(x for x in ref[raw] if x))
            typed = all(ukr[un][raw]) and all(ref[raw])
            tot += typed; ok += typed and h == r
            anchor_rows.append([un, refname, src, ec, "/".join(map(str, h)), "/".join(map(str, r)),
                                "" if not typed else ("match" if h == r else "MISMATCH")])
        print(f"   anchor {un[:24]:24} {ok}/{tot} typed loci match its reference ({src})")

    # ── outputs ────────────────────────────────────────────────────────────────────────
    p_x = out / O_XLSX
    write_table_s3(p_x, order, H["COMB"], amap, lad_comb, flag_rows, flagged, spans)
    with open(out / O_TXT, "w", encoding="utf-8", newline="") as fh:
        hdr = ["Sample Name"] + [f"{ec} {k}" for _, ec in LOCI for k in (1, 2)]
        lines = ["\t".join(hdr)]
        for nm in order:
            lines.append("\t".join([nm] + [str(v) for raw, _ in LOCI for v in sorted(H["COMB"][nm][raw])]))
        fh.write("\n".join(lines) + "\n")
    write_tsv(out / O_DIFF, ["Accession", "Locus", "Raw (bp)", "Barreneche-314 ladder", "Combined ECPGR ladder"], diff)
    write_tsv(out / O_ANCH, ["UKR anchor", "Reference name", "Reference set", "Locus", "Harmonised UKR", "Reference", "Result"], anchor_rows)
    write_tsv(out / O_MAP, ["Locus", "UKR raw (bp)", "ECPGR (bp)", "Calibration span"],
              [[ec, ua, amap[raw][ua], f"{spans[ec][0]}-{spans[ec][1]}"] for raw, ec in LOCI for ua in sorted(amap[raw])])

    got = {
        "n_accessions": len(order),
        "n_anchors_used": len(used),
        "map_sizes": {ec: len(amap[raw]) for raw, ec in LOCI},
        "spans": spans,
        "ladder_sizes_combined": {ec: len(set(lad_comb[ec])) for _, ec in LOCI},
        "ladder_sizes_b314": {ec: len(set(lad_b314[ec])) for _, ec in LOCI},
        "ladder_diff_accessions": len({d[0] for d in diff}),
        "ladder_diff_cells": len(diff),
        "flagged_accessions": len(flagged),
        "flagged_alleles": len(flag_rows),
    }
    checks = {k: (got[k] == v) for k, v in EXPECTED.items()}
    print("counts:", json.dumps(got))
    manifest = {
        "script": SCRIPT, "version": VERSION,
        "run_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "openpyxl": __import__("openpyxl").__version__,
        "inputs": {k: {"file": p.name, "sha256": sha256(p)} for k, p in inputs.items()},
        "outputs": {p.name: sha256(p) for p in [p_x, out / O_TXT, out / O_DIFF, out / O_ANCH, out / O_MAP]},
        "anchors_used": [{"ukr": u, "reference": r, "set": s} for u, r, s in used],
        "counts": got, "expected": EXPECTED, "checks": checks,
    }
    (out / O_MANI).write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"written to {out}")

    if not a.no_check:
        bad = [k for k, ok in checks.items() if not ok]
        if bad:
            for k in bad:
                print(f"CHECK FAILED  {k}: expected {EXPECTED[k]}, got {got[k]}")
            return 1
        print(f"all {len(checks)} regression checks passed (Table S3 v2, 2026-07-18)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
