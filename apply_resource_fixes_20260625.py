"""
apply_resource_fixes_20260625.py
================================
One-time curation script for INS-1666 (Resource update June 2026).

Reads  : data/01_intermediate/resources/resources_2026-05-14.tsv   (baseline)
Reads  : data/01_intermediate/resources/INS_ResourceFixes_20260625.tsv (fixes)
Writes : data/01_intermediate/resources/resources_2026-06-25.tsv   (output)

Every edit is either:
  - driven directly from the Fixes TSV (URL updates), or
  - hardcoded below with a comment referencing the Fixes TSV row
    (removals, consolidations, typo/content fixes).

The output is designed to be fed straight into package_resources.py.
"""

import csv
import os
import sys

# -- Paths ----------------------------------------------------------------
BASE_DIR = os.path.join("data", "01_intermediate", "resources")
BASELINE_FILE = os.path.join(BASE_DIR, "resources_2026-05-14.tsv")
FIXES_FILE = os.path.join(BASE_DIR, "INS_ResourceFixes_20260625.tsv")
OUTPUT_FILE = os.path.join(BASE_DIR, "resources_2026-06-25.tsv")

# -- Column order (must match package_resources.py EXPECTED_COLUMNS) ------
COLUMNS = [
    "type",
    "resource_uuid",
    "resource_source_id",
    "resource_title",
    "resource_short_description",
    "resource_source_url",
    "resource_tool_type",
    "resource_tool_subtype",
    "resource_research_area",
    "resource_research_type",
    "resource_access",
    "resource_doc",
    "resource_poc_email",
    "resource_poc_name",
    "resource_full_description",
]


# =========================================================================
# 1. IDs TO REMOVE
#    Source: Fixes TSV rows where Issue Type = "Remove"
# =========================================================================
REMOVE_IDS = {
    "cbiit_catalog_001",  # Duplicate of r4r_0049
    "cbiit_catalog_011",  # Consolidate into r4r_0073
    "r4r_0026",           # Older duplicate of r4r_0273
    "r4r_0121",           # Consolidate into r4r_0023
    "r4r_0147",           # Consolidate into r4r_0023
    "r4r_0251",           # Project retired; redirects to CCG main page
    "r4r_0263",           # Older duplicate of r4r_0270
    "r4r_0269",           # Older duplicate of r4r_0274
}


# =========================================================================
# 2. CONSOLIDATION OVERRIDES
#    Approved in chat review.  Merges metadata from absorbed records
#    into the surviving record, then the absorbed IDs are removed above.
#
#    r4r_0023  absorbs  r4r_0121 + r4r_0147
#    r4r_0073  absorbs  cbiit_catalog_011
# =========================================================================
CONSOLIDATIONS = {
    # ----- The Cancer Proteome Atlas (r4r_0121 + r4r_0147 → r4r_0023) ----
    "r4r_0023": {
        "resource_title": "The Cancer Proteome Atlas (TCPA)",
        "resource_short_description": (
            "TCPA is a comprehensive resource for accessing, visualizing, "
            "and analyzing functional proteomics data of patient tumor and "
            "cancer cell line samples."
        ),
        "resource_source_url": "https://tcpaportal.org",
        "resource_tool_type": "Analysis Tools;Datasets and Databases",
        "resource_tool_subtype": "Data Visualization",
        "resource_research_area": "Cancer Biology;Cancer Omics;Cancer Treatment",
        "resource_research_type": "Basic;Translational",
        "resource_access": "Open Access",
        "resource_doc": (
            "Center for Cancer Genomics (CCG);"
            "Informatics Technology for Cancer Research (ITCR)"
        ),
        "resource_poc_email": "hliang1@mdanderson.org;tcpa_info@mdanderson.org",
        "resource_poc_name": "Han Liang",
        "resource_full_description": (
            "TCPA is a comprehensive resource for accessing, visualizing, "
            "and analyzing functional proteomics data of patient tumor and "
            "cancer cell line samples. It includes proteomic, genomic, "
            "transcriptomic, and drug screening data primarily derived from "
            "The Cancer Genome Atlas (TCGA) and related initiatives. TCPA "
            "provides interactive modules for data summary, protein "
            "exploration, visualization, analysis, and download."
        ),
    },
    # ----- Personalized Cancer Therapy (cbiit_catalog_011 → r4r_0073) ----
    "r4r_0073": {
        "resource_title": "Personalized Cancer Therapy",
        "resource_short_description": (
            "The Personalized Cancer Therapy website is a knowledge base "
            "for physicians and patients to assess potential therapy options "
            "based on specific tumor biomarkers and genomic alterations, "
            "regardless of disease site."
        ),
        "resource_source_url": "https://pct.mdanderson.org/",
        "resource_tool_type": "Analysis Tools;Datasets and Databases",
        "resource_tool_subtype": "Genomic Analysis;Genomic Datasets",
        "resource_research_area": "Cancer Omics;Cancer Treatment",
        "resource_research_type": "Clinical;Translational",
        "resource_access": "Open Access",
        "resource_doc": "Informatics Technology for Cancer Research (ITCR)",
        "resource_poc_email": "PODS@mdanderson.org",
        "resource_poc_name": "",
        "resource_full_description": (
            "Personalized cancer therapy is a treatment strategy centered "
            "on the ability to predict which patients are more likely to "
            "respond to specific cancer therapies. This approach is founded "
            "upon the idea that tumor biomarkers are associated with patient "
            "prognosis and tumor response to therapy. In addition, patient "
            "genetic factors can be associated with drug metabolism, drug "
            "response and drug toxicity. Personalized tumor molecular "
            "profiles, tumor disease site and other patient characteristics "
            "are then potentially used for determining optimum "
            "individualized therapy options. Tumor biomarkers can be DNA, "
            "RNA, protein and metabolomic profiles that predict therapy "
            "response. However, the most recent approach is the sequencing "
            "of tumor DNA, which can reveal genomic alterations that have "
            "implications for cancer treatment. This Personalized Cancer "
            "Therapy website was specifically developed as a tool for "
            "physicians and patients to assess potential therapy options "
            "based on specific tumor biomarkers."
        ),
    },
}


# =========================================================================
# Helper: read URL fixes straight from the Fixes TSV
# =========================================================================
def load_url_fixes(fixes_path: str) -> dict[str, str]:
    """Return {resource_source_id: new_url} for every row in the Fixes TSV
    where Field == 'website' and New Link is non-empty.

    Rows whose IDs are in REMOVE_IDS are skipped (they'll be deleted).
    Rows whose IDs are in CONSOLIDATIONS are skipped (URL set there).
    """
    # Try utf-8-sig first, fall back to cp1252 (handles Excel artifacts)
    content = None
    for encoding in ("utf-8-sig", "cp1252"):
        try:
            with open(fixes_path, "r", encoding=encoding) as fh:
                content = fh.read()
            break
        except (UnicodeDecodeError, ValueError):
            continue
    if content is None:
        print(f"[FATAL] Cannot decode {fixes_path}")
        sys.exit(1)

    url_fixes: dict[str, str] = {}
    reader = csv.DictReader(content.splitlines(), delimiter="\t")
    for row in reader:
        rid = row.get("Resource ID", "").strip()
        field = row.get("Field", "").strip().lower()
        new_link = row.get("New Link", "").strip()

        if (
            field == "website"
            and new_link
            and rid not in REMOVE_IDS
            and rid not in CONSOLIDATIONS
        ):
            url_fixes[rid] = new_link

    return url_fixes


# =========================================================================
# 3. FIELD-LEVEL FIXES  (content changes, typo corrections)
#    Each block references its Fixes TSV row for traceability.
# =========================================================================
def apply_field_fixes(row: dict) -> list[str]:
    """Apply all non-URL, non-consolidation fixes to a single row.

    Returns a list of human-readable change descriptions for the log.
    """
    rid = row.get("resource_source_id", "")
    changes: list[str] = []

    # ----- r4r_0024 | Typo P1 | "pathothogy" in full_description ---------
    if rid == "r4r_0024":
        old = row["resource_full_description"]
        row["resource_full_description"] = old.replace(
            "pathothogy", "pathology"
        )
        if old != row["resource_full_description"]:
            changes.append("  full_description: 'pathothogy' -> 'pathology'")

    # ----- r4r_0031 | Typo P1 | "dessiminating" in both descriptions -----
    if rid == "r4r_0031":
        for field in ("resource_short_description", "resource_full_description"):
            old = row[field]
            row[field] = old.replace("dessiminating", "disseminating")
            if old != row[field]:
                changes.append(f"  {field}: 'dessiminating' -> 'disseminating'")

    # ----- r4r_0051 | Content P1 | Deduplicate title ---------------------
    if rid == "r4r_0051":
        old = row["resource_title"]
        row["resource_title"] = (
            "Bayesian Phase II Single Arm Clinical Trials "
            "- Time to Event Outcomes Calculator"
        )
        changes.append(f"  title: '{old}' -> '{row['resource_title']}'")

    # ----- r4r_0052 | Content P1 | Deduplicate title ---------------------
    if rid == "r4r_0052":
        old = row["resource_title"]
        row["resource_title"] = (
            "Bayesian Phase II Single Arm Clinical Trials "
            "- Binary Outcomes Calculator"
        )
        changes.append(f"  title: '{old}' -> '{row['resource_title']}'")

    # ----- r4r_0054 | Typo P1 | Pluralize "investigator" -----------------
    if rid == "r4r_0054":
        old = row["resource_short_description"]
        row["resource_short_description"] = old.replace(
            "investigator who use", "investigators who use"
        )
        if old != row["resource_short_description"]:
            changes.append(
                "  short_description: 'investigator who use' "
                "-> 'investigators who use'"
            )

    # ----- r4r_0071 | Typo P1 | "NCI-sponsered" -------------------------
    if rid == "r4r_0071":
        old = row["resource_short_description"]
        row["resource_short_description"] = old.replace(
            "NCI-sponsered", "NCI-sponsored"
        )
        if old != row["resource_short_description"]:
            changes.append(
                "  short_description: 'NCI-sponsered' -> 'NCI-sponsored'"
            )

    # ----- r4r_0080 | Content P2 | Wrong POC email -----------------------
    if rid == "r4r_0080":
        old = row["resource_poc_email"]
        row["resource_poc_email"] = "cronink@mail.nih.gov"
        changes.append(f"  poc_email: '{old}' -> '{row['resource_poc_email']}'")

    # ----- r4r_0082 | Content P2 | Escaped asterisk in HD*Calc -----------
    if rid == "r4r_0082":
        # Short description: HD\*Calc -> HD*Calc
        old_s = row["resource_short_description"]
        row["resource_short_description"] = old_s.replace("\\*", "*")
        if old_s != row["resource_short_description"]:
            changes.append(
                "  short_description: removed backslash-escape from HD*Calc"
            )

        # Full description: HD\<em>Calc...HD\</em>Calc -> HD*Calc...HD*Calc
        old_f = row["resource_full_description"]
        new_f = (
            old_f
            .replace("\\<em>", "*")
            .replace("\\</em>", "*")
            .replace("\\*", "*")
        )
        row["resource_full_description"] = new_f
        if old_f != new_f:
            changes.append(
                "  full_description: removed backslash/em-tag escapes "
                "from HD*Calc"
            )

    # ----- r4r_0085 | Typo P1 | Extra ". " (double period) ---------------
    if rid == "r4r_0085":
        old = row["resource_full_description"]
        row["resource_full_description"] = old.replace(". . ", ". ")
        if old != row["resource_full_description"]:
            changes.append(
                "  full_description: removed extra '. ' (double period/space)"
            )

    # ----- r4r_0098 | Typo P1 | Capitalize title; "resouce" x2 ----------
    if rid == "r4r_0098":
        old_t = row["resource_title"]
        row["resource_title"] = old_t.replace(
            "Growth inhibition Data", "Growth Inhibition Data"
        )
        if old_t != row["resource_title"]:
            changes.append(f"  title: '{old_t}' -> '{row['resource_title']}'")

        old_f = row["resource_full_description"]
        count = old_f.count("resouce")
        row["resource_full_description"] = old_f.replace(
            "resouce", "resource"
        )
        if count:
            changes.append(
                f"  full_description: 'resouce' -> 'resource' ({count}x)"
            )

    # ----- r4r_0130 | Typo P2 | Email domain typo (ni.gov -> nih.gov) ----
    if rid == "r4r_0130":
        old = row["resource_poc_email"]
        row["resource_poc_email"] = old.replace(
            "@mail.ni.gov", "@mail.nih.gov"
        )
        if old != row["resource_poc_email"]:
            changes.append(f"  poc_email: '{old}' -> '{row['resource_poc_email']}'")

    # ----- r4r_0137 | Typo P1 | Convert bare URL to <a href> for dbGaP --
    if rid == "r4r_0137":
        old = row["resource_full_description"]
        row["resource_full_description"] = old.replace(
            "<http://www.ncbi.nlm.nih.gov/gap>",
            '<a href="http://www.ncbi.nlm.nih.gov/gap">'
            "http://www.ncbi.nlm.nih.gov/gap</a>",
        )
        if old != row["resource_full_description"]:
            changes.append(
                "  full_description: converted bare dbGaP URL to <a href>"
            )

    # ----- r4r_0145 | Typo P2 | "nNon-parametric" + backslash -----------
    if rid == "r4r_0145":
        old_s = row["resource_short_description"]
        new_s = old_s.replace("nNon-parametric", "non-parametric")
        new_s = new_s.replace("\\*", "*")
        row["resource_short_description"] = new_s
        if old_s != new_s:
            changes.append(
                "  short_description: 'nNon-parametric' -> 'non-parametric'; "
                "removed backslash from Head*Bang"
            )

        old_f = row["resource_full_description"]
        row["resource_full_description"] = old_f.replace("\\*", "*")
        if old_f != row["resource_full_description"]:
            changes.append(
                "  full_description: removed backslash from Head*Bang"
            )

    # ----- Semicolon-list sort corrections --------------------------------
    # These research_area values are unsorted in the baseline TSV.
    # Hardcode the alphabetical order so package_resources.py report is clean.
    RESEARCH_AREA_SORT = {
        "r4r_0196": "Cancer Biology;Cancer Treatment",
        "r4r_0197": "Bioinformatics;Cancer Biology;Cancer Treatment",
        "r4r_0208": "Cancer Health Disparities;Cancer Statistics;Cancer Treatment",
        "r4r_0229": "Cancer Diagnosis;Cancer Prevention;Screening and Detection",
        "r4r_0230": "Cancer Biology;Cancer Treatment",
        "r4r_0234": "Bioinformatics;Cancer Biology;Cancer Treatment",
        "r4r_0235": "Bioinformatics;Cancer Biology;Cancer Treatment",
        "r4r_0237": "Cancer Biology;Cancer Treatment",
        "r4r_0242": "Cancer Biology;Cancer Omics",
    }
    if rid in RESEARCH_AREA_SORT:
        old = row["resource_research_area"]
        row["resource_research_area"] = RESEARCH_AREA_SORT[rid]
        if old != row["resource_research_area"]:
            changes.append(
                f"  research_area: re-sorted '{old}' -> "
                f"'{row['resource_research_area']}'"
            )

    # ----- r4r_0170 | Content P1 | Set access to Open Access -------------
    if rid == "r4r_0170":
        old = row["resource_access"]
        row["resource_access"] = "Open Access"
        changes.append(f"  access: '{old}' -> 'Open Access'")

    # ----- r4r_0208 | Content P1 | Set access + research_type ------------
    if rid == "r4r_0208":
        old_a = row["resource_access"]
        old_r = row["resource_research_type"]
        row["resource_access"] = "Requires Registration"
        row["resource_research_type"] = "Clinical;Epidemiologic"
        changes.append(f"  access: '{old_a}' -> 'Requires Registration'")
        changes.append(
            f"  research_type: '{old_r}' -> 'Clinical;Epidemiologic'"
        )

    return changes


# =========================================================================
# Main pipeline
# =========================================================================
def main() -> None:
    # -- Verify input files exist -----------------------------------------
    for path in (BASELINE_FILE, FIXES_FILE):
        if not os.path.isfile(path):
            print(f"[FATAL] File not found: {path}")
            sys.exit(1)

    # -- 1. Read baseline -------------------------------------------------
    print(f"Reading baseline: {BASELINE_FILE}")
    with open(BASELINE_FILE, "r", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows: list[dict] = []
        for raw in reader:
            row = {
                k: (v.strip() if v else "")
                for k, v in raw.items()
                if k in COLUMNS
            }
            rows.append(row)
    print(f"  {len(rows)} rows read\n")

    # -- 2. Load URL fixes from Fixes TSV ---------------------------------
    print(f"Reading fixes: {FIXES_FILE}")
    url_fixes = load_url_fixes(FIXES_FILE)
    print(f"  {len(url_fixes)} URL fixes loaded\n")

    # -- 3. Remove rows ---------------------------------------------------
    print(f"--- REMOVALS ({len(REMOVE_IDS)} IDs) ---")
    all_ids = {r["resource_source_id"] for r in rows}
    before = len(rows)
    rows = [r for r in rows if r["resource_source_id"] not in REMOVE_IDS]
    removed_count = before - len(rows)
    for rid in sorted(REMOVE_IDS):
        status = "removed" if rid in all_ids else "NOT FOUND in baseline"
        print(f"  {rid}: {status}")
    print(f"  {removed_count} rows removed ({len(rows)} remaining)\n")

    # -- 4. Apply consolidation overrides ---------------------------------
    print(f"--- CONSOLIDATIONS ({len(CONSOLIDATIONS)} records) ---")
    for row in rows:
        rid = row["resource_source_id"]
        if rid in CONSOLIDATIONS:
            overrides = CONSOLIDATIONS[rid]
            print(f"  {rid}: overwriting {len(overrides)} fields")
            for field, value in overrides.items():
                row[field] = value

    print()

    # -- 5. Apply URL fixes from Fixes TSV --------------------------------
    print(f"--- URL FIXES ({len(url_fixes)} links) ---")
    applied_urls = 0
    for row in rows:
        rid = row["resource_source_id"]
        if rid in url_fixes:
            old = row["resource_source_url"]
            row["resource_source_url"] = url_fixes[rid]
            print(f"  {rid}: {old}")
            print(f"       -> {url_fixes[rid]}")
            applied_urls += 1

    remaining_ids = {r["resource_source_id"] for r in rows}
    unmatched = set(url_fixes.keys()) - remaining_ids
    if unmatched:
        print(f"\n  [WARN] URL fixes not matched to any row:")
        for rid in sorted(unmatched):
            print(f"    {rid}: {url_fixes[rid]}")
    print(f"  {applied_urls} URLs updated\n")

    # -- 6. Apply field-level fixes ---------------------------------------
    print("--- FIELD FIXES ---")
    total_field_fixes = 0
    for row in rows:
        changes = apply_field_fixes(row)
        if changes:
            rid = row["resource_source_id"]
            print(f"  {rid}:")
            for c in changes:
                print(f"  {c}")
            total_field_fixes += len(changes)
    print(f"  {total_field_fixes} field-level changes applied\n")

    # -- 7. Write output --------------------------------------------------
    print("--- OUTPUT ---")
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows)} rows -> {OUTPUT_FILE}")
    print(f"\nDone. Feed this file into package_resources.py.")


if __name__ == "__main__":
    main()
