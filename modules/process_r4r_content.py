"""
process_r4r_content.py
----------------------
One-time conversion script that reads R4R (Resources for Researchers)
markdown files from data/00_input/r4r/, parses their YAML front-matter
and markdown body, and writes a flattened, INS-ready TSV to
data/01_intermediate/r4r_resources_conversion.tsv.

This script is designed to be run once to produce the initial TSV from
the original markdown source files. Future resource additions and edits
should be made directly in the TSV; re-running this script would
overwrite those changes.

A conversion report with value distributions, blank-field audit, and
data quality notes is saved to reports/r4r/r4r_conversion_report.txt.

Source markdown structure
-------------------------
Each .md file contains YAML front-matter delimited by --- and a
markdown body below it. The YAML fields are:

  id              - int, unique resource identifier
  title           - str, resource name
  description     - str, short description (YAML block scalar)
  website         - str, resource URL
  toolTypes       - list of {toolType: "category/subcategory"}
  researchAreas   - list of {researchArea: "area_name"}
  researchTypes   - list of {researchType: "type_name"}
  resourceAccess  - dict with {type: "open" | "register" | "cost"}
  docs            - list of {doc: "nci_division_code"}
  poc             - list of contacts (optional); each may have
                    email and name {firstname, lastname}

Output TSV columns
------------------
  type                       - record type (hardcoded "resource")
  resource_uuid              - deterministic UUID5 from id + title
  resource_source_id         - original numeric id
  resource_title             - cleaned title text
  resource_short_description - short description (markdown -> HTML)
  resource_source_url        - resource website URL
  resource_tool_type         - semicolon-separated parent categories
  resource_tool_subtype      - semicolon-separated subcategories
  resource_research_area     - semicolon-separated research areas
  resource_research_type     - semicolon-separated research types
  resource_access            - access level (open / register / cost)
  resource_doc               - semicolon-separated NCI division codes
  resource_poc_email         - semicolon-separated contact emails
  resource_poc_name          - semicolon-separated contact names
  resource_full_description  - full narrative text (markdown -> HTML)


Processing steps
----------------
  1. Parse YAML front-matter and markdown body from each .md file
  2. Flatten nested YAML structures into semicolon-separated strings
  3. Split toolTypes at "/" into parent type and subtype columns
  4. Convert markdown formatting (bold, italic, links, lists) to HTML
  5. Sanitise text for TSV safety (collapse whitespace, strip)
  6. Generate deterministic UUID5 for each resource
  7. Write TSV and conversion report

Post-generation curation
------------------------
After initial TSV generation, AI-recommended values (via GitHub Copilot)
were applied to resources that had blank fields for short_description,
tool_type, or research_area.  Curated values live in the CURATIONS dict
below.  See reports/r4r/r4r_curation_changelog.txt for the full change log.
"""

import os
import glob
import csv
import re
import uuid
from collections import Counter 
from datetime import datetime
import yaml
import unicodedata

# -- UUID namespace (matches convention in package_output_data.py) --------
R4R_UUID_NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')

# -- Paths ----------------------------------------------------------------
INPUT_DIR   = os.path.join("data", "00_input", "r4r")
OUTPUT_TSV  = os.path.join("data", "01_intermediate", "r4r_resources_conversion.tsv")
REPORT_PATH = os.path.join("reports", "r4r", "r4r_conversion_report.txt")

# -- Output column order -------------------------------------------------
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

# -- Human-friendly label overrides ---------------------------------------
# The default transform is:  value.replace("_", " ").title()
# Entries here override that default for values that need special casing,
# punctuation, or wording.  Add new overrides as needed.
LABEL_OVERRIDES: dict[str, str] = {
    # tool types
    "datasets_databases":       "Datasets and Databases",
    "networks_consortiums":     "Networks and Consortiums",
    # tool subtypes
    "guidelines_protocols":     "Guidelines and Protocols",
    # research areas
    "causes_of_cancer":         "Causes of Cancer",
    "screening_detection":      "Screening and Detection",
    "cancer_omics":             "Cancer Omics",
}


def humanize_label(raw_value: str) -> str:
    """Convert a snake_case metadata value to a human-friendly label.

    Checks LABEL_OVERRIDES first; falls back to replacing underscores
    with spaces and applying title case.
    """
    if not raw_value:
        return raw_value
    if raw_value in LABEL_OVERRIDES:
        return LABEL_OVERRIDES[raw_value]
    return raw_value.replace("_", " ").title()


def humanize_semicolon_list(raw_field: str) -> str:
    """Apply humanize_label to each value in a semicolon-separated string."""
    if not raw_field:
        return raw_field
    return ";".join(humanize_label(v.strip()) for v in raw_field.split(";"))


# -- Access-level mapping ------------------------------------------------
ACCESS_LABELS: dict[str, str] = {
    "open":     "Open Access",
    "register": "Requires Registration",
    "cost":     "Requires Payment",
}


def humanize_access(raw_value: str) -> str:
    """Map a raw access value to its human-friendly label."""
    return ACCESS_LABELS.get(raw_value.strip(), raw_value)


def dedup_semicolon_list(raw_field: str) -> str:
    """Remove duplicate values from a semicolon-separated string, preserving order.

    Empty tokens (e.g. from trailing/double semicolons) are also filtered out.
    """
    if not raw_field:
        return raw_field
    seen: set[str] = set()
    result: list[str] = []
    for v in raw_field.split(";"):
        v = v.strip()
        if v and v not in seen:
            seen.add(v)
            result.append(v)
    return ";".join(result)


def sort_semicolon_list(raw_field: str) -> str:
    """Sort values in a semicolon-separated string alphabetically (case-insensitive).

    Empty/blank tokens are filtered out.
    """
    if not raw_field:
        return raw_field
    values = [v.strip() for v in raw_field.split(";") if v.strip()]
    return ";".join(sorted(values, key=str.casefold))


# -- NCI division / office code mapping ----------------------------------
DOC_LABELS: dict[str, str] = {
    "itcr":   "Informatics Technology for Cancer Research (ITCR)",
    "dctd":   "Division of Cancer Treatment & Diagnosis (DCTD)",
    "ccg":    "Center for Cancer Genomics (CCG)",
    "dccps":  "Division of Cancer Control & Population Sciences (DCCPS)",
    "dceg":   "Division of Cancer Epidemiology & Genetics (DCEG)",
    "ocg":    "Office of Cancer Genomics (OCG)",
    "dcb":    "Division of Cancer Biology (DCB)",
    "cbiit":  "Center for Biomedical Informatics & Information Technology (CBIIT)",
    "cssi":   "Center for Strategic Scientific Initiatives (CSSI)",
    "ocpl":   "Office of Communications & Public Liaison (OCPL)",
    "occpr":  "Office of Cancer Clinical Proteomics Research (OCCPR)",
    "dcp":    "Division of Cancer Prevention (DCP)",
    "oer":    "Office of Extramural Research (OER)",
    "ccr":    "Center for Cancer Research (CCR)",
}


def humanize_doc(raw_value: str) -> str:
    """Map a raw NCI division/office code to its full name."""
    return DOC_LABELS.get(raw_value.strip(), raw_value)


def humanize_doc_list(raw_field: str) -> str:
    """Apply humanize_doc to each value in a semicolon-separated string."""
    if not raw_field:
        return raw_field
    return ";".join(humanize_doc(v.strip()) for v in raw_field.split(";"))


# -- Special character replacement ---------------------------------------
# Mirrors the approach in package_output_data.replace_defined_characters
# but applied directly to R4R string fields so the TSV is clean without
# needing to route through that module.

_CHAR_REPLACEMENTS: dict[int, str] = {
    # Punctuation / whitespace
    0x2028: ' ',       # Line Separator -> space
    0x2029: ' ',       # Paragraph Separator -> space
    0x00A0: ' ',       # Non-breaking space -> space
    0xFEFF: '',        # BOM -> remove
    0x2013: '-',       # En dash -> hyphen
    0x2014: '-',       # Em dash -> hyphen
    0x2010: '-',       # Unicode hyphen -> ASCII hyphen
    0x201C: "'",       # Left double quote -> single quote
    0x201D: "'",       # Right double quote -> single quote
    0x2018: "'",       # Left single quote -> single quote
    0x2019: "'",       # Right single quote -> single quote
    # Greek letters -> HTML entities for proper page rendering
    0x03B1: '&alpha;',   # alpha
    0x03B2: '&beta;',    # beta
    0x03B3: '&gamma;',   # gamma
    0x03B4: '&delta;',   # delta (lower)
    0x03B5: '&epsilon;', # epsilon (lower)
    0x03F5: '&epsilon;', # epsilon variant
    0x03BA: '&kappa;',   # kappa
    0x03BB: '&lambda;',  # lambda
    0x03BC: '&mu;',      # mu
    0x0394: '&Delta;',   # Delta (upper)
    0x03A3: '&Sigma;',   # Sigma (upper)
    # Superscripts -> <sup> HTML tags for proper page rendering
    0x2070: '<sup>0</sup>', 0x00B9: '<sup>1</sup>',
    0x00B2: '<sup>2</sup>', 0x00B3: '<sup>3</sup>',
    0x2074: '<sup>4</sup>', 0x2075: '<sup>5</sup>',
    0x2076: '<sup>6</sup>', 0x2077: '<sup>7</sup>',
    0x2078: '<sup>8</sup>', 0x2079: '<sup>9</sup>',
    # Subscripts -> <sub> HTML tags for proper page rendering
    0x2080: '<sub>0</sub>', 0x2081: '<sub>1</sub>',
    0x2082: '<sub>2</sub>', 0x2083: '<sub>3</sub>',
    0x2084: '<sub>4</sub>', 0x2085: '<sub>5</sub>',
    0x2086: '<sub>6</sub>', 0x2087: '<sub>7</sub>',
    0x2088: '<sub>8</sub>', 0x2089: '<sub>9</sub>',
    # Symbols -> HTML entities for proper page rendering
    0x00A9: '&copy;',     # Copyright
    0x00AE: '&reg;',      # Registered trademark
    0x2122: '&trade;',    # Trademark
    0x00D7: '&times;',    # Multiplication sign
    0x2264: '&le;',       # Less than or equal
    0x2265: '&ge;',       # Greater than or equal
    # Ligatures
    0xFB00: 'ff',      # ff ligature
    0xFB01: 'fi',      # fi ligature
    0xFB02: 'fl',      # fl ligature
    0xFB03: 'ffi',     # ffi ligature
    0xFB04: 'ffl',     # ffl ligature
}

_TRANSLATION_TABLE = str.maketrans(_CHAR_REPLACEMENTS)


def replace_special_characters(text: str) -> str:
    """Replace non-standard characters with ASCII-safe equivalents.

    Uses a translation table modeled after package_output_data's
    replace_defined_characters.  Any remaining non-ASCII characters
    (after translation) are dropped via NFD normalization + ASCII
    encoding, which also strips accents from Latin characters.
    """
    if not text:
        return text

    # Apply explicit replacements
    text = text.translate(_TRANSLATION_TABLE)

    # Drop any remaining non-ASCII (accents, stray Unicode, etc.)
    text = (unicodedata.normalize('NFD', text)
            .encode('ascii', 'ignore')
            .decode('utf-8'))

    return text


def clean_for_tsv(text: str) -> str:
    """Collapse whitespace so a string is safe inside a single TSV cell.

    Used for short metadata fields (e.g. resource_title) where line
    breaks are not meaningful.
    """
    text = text.replace("\r\n", " ").replace("\r", " ").replace("\n", " ")
    text = text.replace("\t", " ")
    text = re.sub(r" {2,}", " ", text)
    return text.strip()


def markdown_to_html(text: str) -> str:
    """Convert lightweight markdown formatting to inline HTML.

    Handles the patterns actually found in the R4R corpus:
      - Markdown links   [text](url)   -> <a href="url">text</a>
      - Bold             **text**      -> <strong>text</strong>
      - Italic           *text*        -> <em>text</em>  (after bold)
      - Unordered lists  * item / - item on their own lines
      - Newlines / paragraph breaks    -> <br>
      - Tabs, multi-spaces             -> collapsed
    Existing HTML tags (<sub>, <br>, etc.) are preserved as-is.
    """
    # Normalise line endings
    text = text.replace("\r\n", "\n").replace("\r", "\n")

    # -- Markdown links -> <a> --------------------------------------------
    text = re.sub(
        r"\[([^\]]+)\]\(([^)]+)\)",
        r'<a href="\2">\1</a>',
        text,
    )

    # -- Bold **text** -> <strong> (before italic) ------------------------
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)

    # -- Italic *text* -> <em> --------------------------------------------
    # Avoid matching list markers (* item) at the start of a line
    text = re.sub(r"(?<!\*)\*(?!\s)(.+?)(?<!\s)\*(?!\*)", r"<em>\1</em>", text)

    # -- Unordered list items  (* item  or  - item) ----------------------
    # Convert markdown list markers at the beginning of a line to a
    # hyphen so they survive the whitespace collapse below.
    text = re.sub(r"(?m)^[*\-+]\s+", "- ", text)

    # -- Collapse whitespace for TSV safety ------------------------------
    # Replace newlines with <br>, tabs with space, multi-spaces with one.
    text = text.replace("\t", " ")
    text = text.replace("\n", "<br>")
    text = re.sub(r"(<br>){2,}", "<br>", text)   # collapse consecutive <br>
    text = re.sub(r" {2,}", " ", text)
    text = text.strip()
    text = re.sub(r"(<br>\s*)+$", "", text)      # strip trailing <br>

    return text


def parse_markdown_file(filepath: str) -> dict:
    """Parse a single R4R markdown file into a flat dict of resource fields."""
    with open(filepath, "r", encoding="utf-8") as fh:
        raw = fh.read()

    # -- Split YAML front-matter from markdown body ----------------------
    # Files are structured as:  ---\n<yaml>\n---\n<body>
    parts = raw.split("---", 2)
    if len(parts) < 3:
        raise ValueError(f"Could not parse front-matter in {filepath}")

    yaml_text = parts[1]
    body_text = parts[2].strip()

    meta = yaml.safe_load(yaml_text)

    # -- Flatten list-of-dict fields into semicolon-separated strings ----
    # toolTypes are split at "/" into parent type and subtype columns
    raw_tool_types = [item["toolType"] for item in (meta.get("toolTypes") or [])]
    tool_types     = ";".join(t.split("/")[0] for t in raw_tool_types)
    tool_subtypes = ";".join(t.split("/", 1)[1] if "/" in t else "" for t in raw_tool_types)

    research_areas = ";".join(
        item["researchArea"] for item in (meta.get("researchAreas") or [])
    )
    research_types = ";".join(
        item["researchType"].strip()
        for item in (meta.get("researchTypes") or [])
    )
    docs = ";".join(
        item["doc"] for item in (meta.get("docs") or [])
    )

    # -- resourceAccess - single-value dict {"type": "open"} --------------
    resource_access = ""
    ra = meta.get("resourceAccess")
    if isinstance(ra, dict):
        resource_access = ra.get("type", "")
    elif isinstance(ra, str):
        resource_access = ra
    resource_access = humanize_access(resource_access)

    # -- poc (point of contact) - optional list of contacts ----------------
    poc_emails = []
    poc_names  = []
    for contact in (meta.get("poc") or []):
        if isinstance(contact, dict):
            email = contact.get("email", "")
            if email:
                poc_emails.append(email)
            name_info = contact.get("name")
            if isinstance(name_info, dict):
                first = (name_info.get("firstname") or "").strip()
                last  = (name_info.get("lastname")  or "").strip()
                full  = f"{first} {last}".strip()
                if full:
                    poc_names.append(full)

    # -- Generate deterministic UUID5 from id + title ---------------------
    r4r_id = str(meta.get("id", ""))
    r4r_title = (meta.get("title") or "").strip()
    resource_uuid = str(uuid.uuid5(
        R4R_UUID_NAMESPACE, "||".join([r4r_id, r4r_title])
    ))

    return {
        "type":                      "resource",
        "resource_uuid":             resource_uuid,
        "resource_source_id":        r4r_id,
        "resource_title":            clean_for_tsv(r4r_title),
        "resource_short_description":      markdown_to_html((meta.get("description") or "")) or "No description preview.",
        "resource_source_url":       (meta.get("website") or "").strip(),
        "resource_tool_type":        sort_semicolon_list(dedup_semicolon_list(humanize_semicolon_list(tool_types))) or "No Tool Type provided",
        "resource_tool_subtype":     sort_semicolon_list(dedup_semicolon_list(humanize_semicolon_list(tool_subtypes))),
        "resource_research_area":    sort_semicolon_list(dedup_semicolon_list(humanize_semicolon_list(research_areas))) or "No Research Area provided",
        "resource_research_type":    sort_semicolon_list(dedup_semicolon_list(humanize_semicolon_list(research_types))),
        "resource_access":           resource_access,
        "resource_doc":              sort_semicolon_list(dedup_semicolon_list(humanize_doc_list(docs))),
        "resource_poc_email":        sort_semicolon_list(dedup_semicolon_list(";".join(poc_emails))),
        "resource_poc_name":         sort_semicolon_list(dedup_semicolon_list(";".join(poc_names))),
        "resource_full_description": markdown_to_html(body_text),
    }


def process_all_r4r_files() -> list[dict]:
    """Find all .md files in INPUT_DIR, parse each, and return rows sorted by id."""
    pattern = os.path.join(INPUT_DIR, "*.md")
    files = sorted(glob.glob(pattern))

    if not files:
        print(f"[WARN] No .md files found in {INPUT_DIR}")
        return []

    rows = []
    for fp in files:
        try:
            row = parse_markdown_file(fp)
            rows.append(row)
        except Exception as exc:
            print(f"[ERROR] Failed to parse {fp}: {exc}")

    # Sort by numeric resource_source_id
    rows.sort(key=lambda r: int(r["resource_source_id"]))
    return rows


def sanitize_rows(rows: list[dict]) -> None:
    """Replace non-standard characters across all string fields in-place.

    Should be called after all other processing (parsing, humanization,
    curation) and before writing the TSV and generating the report.
    """
    # Columns where HTML tags are allowed (description fields only)
    _DESCRIPTION_COLS = {"resource_short_description", "resource_full_description"}

    for row in rows:
        for key, val in row.items():
            if isinstance(val, str):
                row[key] = replace_special_characters(val)

        # Normalize "CTD^2" and "CTD2" to "CTD<sup>2</sup>" in description
        # columns only.  Title columns use plain "CTD2" (set via CURATIONS).
        for col in _DESCRIPTION_COLS:
            if col in row and isinstance(row[col], str):
                row[col] = re.sub(r'CTD\^2', r'CTD<sup>2</sup>', row[col])
                row[col] = re.sub(r'CTD2(?!</)', r'CTD<sup>2</sup>', row[col])


def write_tsv(rows: list[dict], output_path: str) -> None:
    """Write resource rows to a tab-separated file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[OK] Wrote {len(rows)} rows -> {output_path}")


def generate_report(rows: list[dict]) -> str:
    """Build a plain-text conversion report with data quality summary."""
    lines: list[str] = []
    w = lines.append  # shorthand

    w("=" * 70)
    w("  R4R Markdown -> TSV Conversion Report")
    w(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    w("=" * 70)

    # -- Overview --------------------------------------------------------
    int_ids = sorted(int(r["resource_source_id"]) for r in rows)
    w("")
    w("OVERVIEW")
    w(f"  Total resources processed : {len(rows)}")
    w(f"  ID range                  : {int_ids[0]} - {int_ids[-1]}")
    expected = set(range(int_ids[0], int_ids[-1] + 1))
    missing_ids = sorted(expected - set(int_ids))
    w(f"  Missing IDs in range      : {missing_ids if missing_ids else 'None'}")
    w(f"  Input directory            : {INPUT_DIR}")
    w(f"  Output TSV                 : {OUTPUT_TSV}")
    w(f"  Columns ({len(COLUMNS)})              : {', '.join(COLUMNS)}")

    # -- Blank / empty field audit ---------------------------------------
    REQUIRED = ["resource_source_id", "resource_title", "resource_short_description",
                "resource_source_url", "resource_tool_type",
                "resource_tool_subtype", "resource_research_area",
                "resource_research_type", "resource_access", "resource_doc"]
    w("")
    w("BLANK FIELD AUDIT")
    any_blanks = False
    for col in REQUIRED:
        blanks = [str(r["resource_source_id"]) for r in rows if not str(r.get(col, "")).strip()]
        if blanks:
            any_blanks = True
            preview = ", ".join(blanks[:15])
            suffix = f"  ... and {len(blanks) - 15} more" if len(blanks) > 15 else ""
            w(f"  {col:<30s}  {len(blanks):>3} blank  -  ids: {preview}{suffix}")
    if not any_blanks:
        w("  All required fields populated.")

    # -- Value distributions ---------------------------------------------
    w("")
    w("VALUE DISTRIBUTIONS")
    for field in ("resource_tool_type", "resource_tool_subtype",
                  "resource_research_area", "resource_research_type",
                  "resource_access", "resource_doc"):
        vals: Counter = Counter()
        for r in rows:
            for v in r[field].split(";"):
                v = v.strip()
                if v:
                    vals[v] += 1
        w(f"\n  {field}  ({len(vals)} unique values):")
        for v, c in vals.most_common():
            w(f"      {c:>4}  {v}")

    # -- POC coverage ----------------------------------------------------
    poc_with_email    = sum(1 for r in rows if r["resource_poc_email"].strip())
    poc_with_name     = sum(1 for r in rows if r["resource_poc_name"].strip())
    poc_email_no_name = sum(1 for r in rows
                           if r["resource_poc_email"].strip() and not r["resource_poc_name"].strip())
    w("")
    w("POINT-OF-CONTACT (POC) COVERAGE")
    w(f"  With email          : {poc_with_email}/{len(rows)}")
    w(f"  With name           : {poc_with_name}/{len(rows)}")
    w(f"  Email but no name   : {poc_email_no_name}")

    # -- Description length stats ----------------------------------------
    desc_lens = [len(r["resource_full_description"]) for r in rows]
    empty_desc = [str(r["resource_source_id"]) for r in rows if not r["resource_full_description"].strip()]
    w("")
    w("FULL DESCRIPTION STATS")
    w(f"  Min length : {min(desc_lens)} chars")
    w(f"  Max length : {max(desc_lens)} chars")
    w(f"  Avg length : {sum(desc_lens) // len(desc_lens)} chars")
    w(f"  Empty      : {', '.join(empty_desc) if empty_desc else 'None'}")

    # -- Preview ---------------------------------------------------------
    w("")
    w("PREVIEW (first 10 rows)")
    w(f"  {'id':>4}  {'title':<55s}  {'tool_types'}")
    w(f"  {'--':>4}  {'-' * 55}  {'-' * 40}")
    for r in rows[:10]:
        w(f"  {r['resource_source_id']:>4}  {r['resource_title'][:55]:<55s}  {r['resource_tool_type']}")

    w("")
    return "\n".join(lines)


def write_report(report_text: str, report_path: str) -> None:
    """Save the conversion report to a text file."""
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report_text + "\n")
    print(f"[OK] Report saved -> {report_path}")


# -- Post-generation curation --------------------------------------------
# Curated and AI-recommended values (GitHub Copilot) for resources whose 
# original markdown files needed edits or had blank fields. 
# Fixes are applied in-memory before writing the TSV.
# See reports/r4r/r4r_curation_changelog.txt for the full audit trail.

CURATIONS: dict[str, dict[str, str]] = {
    # -- resource_title typo fixes (4 resources) ----------------------------------------
    "28": {"resource_title": "Tumor Heterogeneity Research Interactive Visualization Environment (THRIVE)"}, # Heterogenity -> Heterogeneity
    "47": {"resource_title": "Pathological Complete Response (pCR) Trial-Level Surrogate Analysis Software"}, # Added missing space
    "51": {"resource_title": "Bayesian Phase II Single Arm Clinical Trials"}, # Baysian -> Bayesian
    "52": {"resource_title": "Bayesian Phase II Single Arm Clinical Trials"}, # Baysian -> Bayesian
    "126": {"resource_title": "Cancer Target Discovery and Development (CTD2) Dashboard"}, # CTD^2 -> CTD2 (plain text, no HTML in titles)
    "135": {"resource_title": "Cancer Target Discovery and Development (CTD2) Data Portal"}, # CTD^2 -> CTD2 (plain text, no HTML in titles)
    # -- resource_short_description (1 resource) -------------------------
    "44": {
        "resource_short_description":
            "A Microsoft Excel-based program for dose level assignment "
            "in accelerated titration designs for phase I clinical trials.",
    },
    # -- resource_tool_type (6 resources) --------------------------------
    "171": {"resource_tool_type": "Analysis Tools",
            "resource_tool_subtype": "Statistical Software"},
    "172": {"resource_tool_type": "Analysis Tools",
            "resource_tool_subtype": "Statistical Software"},
    "173": {"resource_tool_type": "Community Research Tools",
            "resource_tool_subtype": "Questionnaire"},
    "174": {"resource_tool_type": "Analysis Tools",
            "resource_tool_subtype": "Modeling"},
    "211": {"resource_tool_type": "Clinical Research Tools",
            "resource_tool_subtype": "Guidelines and Protocols",
            "resource_research_area": "Cancer Treatment"},
    "215": {"resource_tool_type": "Community Research Tools",
            "resource_research_area": "Cancer Biology;Cancer Treatment"},
    # -- resource_research_area (remaining 30 resources) -----------------
    "62":  {"resource_research_area": "Bioinformatics"},
    "63":  {"resource_research_area": "Bioinformatics"},
    "64":  {"resource_research_area": "Bioinformatics"},
    "182": {"resource_research_area": "Cancer Prevention"},
    "193": {"resource_research_area": "Cancer Biology;Cancer Omics"},
    "194": {"resource_research_area": "Cancer Biology;Cancer Treatment"},
    "195": {"resource_research_area": "Cancer Treatment"},
    "196": {"resource_research_area": "Cancer Treatment;Cancer Biology"},
    "197": {"resource_research_area": "Cancer Biology;Cancer Treatment;Bioinformatics"},
    "201": {"resource_research_area": "Cancer Treatment"},
    "202": {"resource_research_area": "Cancer Biology"},
    "203": {"resource_research_area": "Cancer Treatment"},
    "204": {"resource_research_area": "Cancer Treatment"},
    "205": {"resource_research_area": "Cancer Treatment"},
    "207": {"resource_research_area": "Cancer Treatment"},
    "208": {"resource_research_area": "Cancer Statistics;Cancer Treatment;Cancer Health Disparities"},
    "209": {"resource_research_area": "Cancer Treatment"},
    "210": {"resource_research_area": "Cancer Treatment"},
    "229": {"resource_research_area": "Cancer Prevention;Cancer Diagnosis;Screening and Detection"},
    "230": {"resource_research_area": "Cancer Treatment;Cancer Biology"},
    "233": {"resource_research_area": "Cancer Biology;Cancer Omics"},
    "234": {"resource_research_area": "Cancer Biology;Cancer Treatment;Bioinformatics"},
    "235": {"resource_research_area": "Cancer Biology;Cancer Treatment;Bioinformatics"},
    "236": {"resource_research_area": "Cancer Treatment"},
    "237": {"resource_research_area": "Cancer Treatment;Cancer Biology"},
    "238": {"resource_research_area": "Cancer Biology;Cancer Omics;Causes of Cancer"},
    "239": {"resource_research_area": "Cancer Biology;Cancer Omics"},
    "241": {"resource_research_area": "Cancer Biology;Cancer Treatment"},
    "242": {"resource_research_area": "Cancer Omics;Cancer Biology"},
    "263": {"resource_research_area": "Causes of Cancer;Cancer Omics"},
}

CURATION_LOG_PATH = os.path.join("reports", "r4r", "r4r_curation_changelog.txt")


def apply_curations(rows: list[dict]) -> int:
    """Apply AI-curated values to in-memory rows and write a changelog.

    Returns the number of field-level changes applied.
    """
    changelog: list[str] = []
    changelog.append("=" * 70)
    changelog.append("  R4R Post-Conversion Curation Changelog")
    changelog.append(f"  Applied: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    changelog.append("  Method: AI-recommended values (GitHub Copilot)")
    changelog.append("  Source: Resource titles, full descriptions, and contextual analysis")
    changelog.append("=" * 70)
    changelog.append("")
    changelog.append("All values below were recommended by an LLM based on each resource's")
    changelog.append("title, full description, and related metadata. These replace placeholder")
    changelog.append("values that were inserted during the initial markdown-to-TSV conversion")
    changelog.append("for resources where the original markdown files had no value for the field.")
    changelog.append("")
    changelog.append("Fields curated:")
    changelog.append("  - resource_short_description  (1 resource)")
    changelog.append("  - resource_tool_type          (6 resources, some with tool_subtype)")
    changelog.append("  - resource_research_area      (32 resources)")
    changelog.append("")
    changelog.append("-" * 70)

    changes = 0
    for row in rows:
        rid = row["resource_source_id"]
        if rid in CURATIONS:
            for field, new_val in CURATIONS[rid].items():
                old_val = row[field]
                row[field] = new_val
                changes += 1
                changelog.append("")
                changelog.append(f"  Resource ID : {rid}")
                changelog.append(f"  Title       : {row['resource_title']}")
                changelog.append(f"  Field       : {field}")
                changelog.append(f"  Old value   : {old_val}")
                changelog.append(f"  New value   : {new_val}")
                changelog.append(f"  {'-' * 60}")

    changelog.append("")
    changelog.append(f"Total changes applied: {changes}")

    os.makedirs(os.path.dirname(CURATION_LOG_PATH), exist_ok=True)
    with open(CURATION_LOG_PATH, "w", encoding="utf-8") as fh:
        fh.write("\n".join(changelog) + "\n")
    print(f"[OK] Applied {changes} curated values -> {CURATION_LOG_PATH}")
    return changes


# -- Main -----------------------------------------------------------------
if __name__ == "__main__":
    rows = process_all_r4r_files()

    # Apply AI-curated values for blank fields
    apply_curations(rows)

    # Replace non-standard characters with ASCII / HTML entities
    sanitize_rows(rows)

    # Write TSV
    write_tsv(rows, OUTPUT_TSV)

    # Generate, print, and save report
    report = generate_report(rows)
    print(report)
    write_report(report, REPORT_PATH)
