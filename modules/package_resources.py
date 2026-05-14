"""
package_resources.py
--------------------
Reads the most recent curated resources_{date}.tsv from
data/01_intermediate/resources/, stamps deterministic UUID5s,
validates data quality, and writes:

  - data/02_output/resources/resources_{date}_processed.tsv
  - reports/resources/resources_{date}_issue_report.txt

Designed to be run after each curation cycle.  The intermediate TSV
may be edited in Excel, so this script handles common Excel artifacts
(BOM, encoding fallback, whitespace, quote escaping).

Exit codes:
  0 - Success (output written, issues may be logged but non-critical)
  1 - Critical error (duplicate UUIDs, missing expected columns)
"""

import os
import sys
import csv
import re
import uuid
import glob
from collections import Counter
from datetime import datetime

# -- Shared UUID namespace  -----------------------------------------------
NAMESPACE = uuid.UUID("12345678-1234-5678-1234-567812345678")

# -- Paths ----------------------------------------------------------------
INPUT_DIR = os.path.join("data", "01_intermediate", "resources")
OUTPUT_DIR = os.path.join("data", "02_output", "resources")
REPORT_DIR = os.path.join("reports", "resources")

# -- Expected columns  ----------------------------------------------------
EXPECTED_COLUMNS = [
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

# -- Required fields (blank = flagged issue) ------------------------------
REQUIRED_FIELDS = [
    "resource_title",
    "resource_source_url",
    "resource_tool_type",
    "resource_research_area",
]

# -- Common non-ASCII characters and their recommended replacements -------
NON_ASCII_SUGGESTIONS: dict[str, str] = {
    "\ufffd": "(encoding error - retype the character manually)",
    "\u2018": "' (left single quote -> apostrophe)",
    "\u2019": "' (right single quote -> apostrophe)",
    "\u201c": "' (left double quote -> apostrophe or remove)",
    "\u201d": "' (right double quote -> apostrophe or remove)",
    "\u2013": "- (en dash -> hyphen)",
    "\u2014": "- (em dash -> hyphen)",
    "\u00a0": "(non-breaking space -> regular space)",
    "\u00b2": "write as <sup>2</sup> in descriptions",
    "\u00b3": "write as <sup>3</sup> in descriptions",
    "\u03b1": "write as &alpha; in descriptions",
    "\u03b2": "write as &beta; in descriptions",
    "\u03bc": "write as &mu; in descriptions",
}

# -- ID format pattern ----------------------------------------------------
ID_PATTERN = re.compile(r"^[a-z][a-z0-9]*(_[a-z0-9]+)+$")


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------


def find_latest_input(input_dir: str = INPUT_DIR) -> str | None:
    """Find the most recent resources_{date}.tsv in the input directory.

    Matches files named resources_YYYY-MM-DD.tsv and returns the one
    with the latest date in the filename.
    """
    pattern = os.path.join(input_dir, "resources_*.tsv")
    candidates = glob.glob(pattern)

    # Filter to those matching the expected date pattern
    dated = []
    for path in candidates:
        basename = os.path.basename(path)
        match = re.match(r"resources_(\d{4}-\d{2}-\d{2})\.tsv$", basename)
        if match:
            dated.append((match.group(1), path))

    if not dated:
        return None

    # Sort by date string (YYYY-MM-DD sorts lexicographically)
    dated.sort(key=lambda x: x[0], reverse=True)
    return dated[0][1]


def extract_date_from_path(filepath: str) -> str:
    """Extract the YYYY-MM-DD date string from a resources_{date}.tsv filename."""
    basename = os.path.basename(filepath)
    match = re.search(r"(\d{4}-\d{2}-\d{2})", basename)
    if match:
        return match.group(1)
    raise ValueError(f"Cannot extract date from filename: {basename}")


def read_resources(filepath: str) -> tuple[list[dict], list[str]]:
    """Read a resources TSV, handling Excel encoding quirks.

    Returns (rows, extra_columns) where extra_columns lists any
    column names not in EXPECTED_COLUMNS.

    Raises SystemExit if any expected columns are completely missing.
    """
    # Try utf-8-sig first (handles BOM), fall back to cp1252
    content = None
    for encoding in ("utf-8-sig", "cp1252"):
        try:
            with open(filepath, "r", encoding=encoding) as fh:
                content = fh.read()
            break
        except (UnicodeDecodeError, ValueError):
            continue

    if content is None:
        print(f"[FATAL] Expected utf-8 or cp1252 encoding for {filepath}.")
        sys.exit(1)

    # Parse TSV
    reader = csv.DictReader(content.splitlines(), delimiter="\t")
    input_columns = list(reader.fieldnames or [])

    # Check for missing expected columns
    missing_cols = [c for c in EXPECTED_COLUMNS if c not in input_columns]
    if missing_cols:
        print(f"[FATAL] Missing expected columns in {filepath}:")
        for col in missing_cols:
            print(f"  - {col}")
        sys.exit(1)

    # Identify extra columns
    extra_cols = [c for c in input_columns if c not in EXPECTED_COLUMNS]

    # Read rows, stripping whitespace from all values
    rows: list[dict] = []
    for raw_row in reader:
        row = {k: (v.strip() if v else "") for k, v in raw_row.items()
               if k in EXPECTED_COLUMNS}
        rows.append(row)

    return rows, extra_cols


def generate_uuids(rows: list[dict]) -> None:
    """Generate deterministic UUID5 for each row based on resource_source_id.

    Overwrites the resource_uuid column in-place.
    """
    for row in rows:
        source_id = row.get("resource_source_id", "").strip()
        if source_id:
            row["resource_uuid"] = str(uuid.uuid5(NAMESPACE, source_id))
        else:
            row["resource_uuid"] = ""


def fix_type_column(rows: list[dict]) -> list[str]:
    """Ensure the type column is 'resource' for every row.

    Returns a list of issue strings for rows where it was incorrect.
    """
    issues: list[str] = []
    for row in rows:
        val = row.get("type", "").strip()
        if not val:
            row["type"] = "resource"
        elif val != "resource":
            issues.append(
                f"Row '{row.get('resource_source_id', '?')}': "
                f"type was '{val}', auto-corrected to 'resource'"
            )
            row["type"] = "resource"
    return issues


def validate(rows: list[dict], extra_columns: list[str]) -> list[str]:
    """Run all validation checks and return a list of issue strings.

    Does NOT include duplicate UUID check (that's handled separately
    as a critical error).
    """
    issues: list[str] = []

    # -- Extra columns ---------------------------------------------------
    if extra_columns:
        issues.append("EXTRA COLUMNS (not in expected schema, ignored in output):")
        for col in extra_columns:
            issues.append(f"  - {col}")
        issues.append("")

    # -- Duplicate resource_source_id ------------------------------------
    id_counts: dict[str, int] = {}
    for row in rows:
        rid = row.get("resource_source_id", "")
        id_counts[rid] = id_counts.get(rid, 0) + 1
    dupes = {k: v for k, v in id_counts.items() if v > 1}
    if dupes:
        issues.append("DUPLICATE resource_source_id VALUES:")
        for rid, count in sorted(dupes.items()):
            issues.append(f"  '{rid}' appears {count} times")
        issues.append("")

    # -- ID format check -------------------------------------------------
    bad_ids: list[str] = []
    for row in rows:
        rid = row.get("resource_source_id", "").strip()
        if rid and not ID_PATTERN.match(rid):
            bad_ids.append(rid)
    if bad_ids:
        issues.append("IDs NOT MATCHING EXPECTED FORMAT (prefix_NNNN):")
        for rid in bad_ids[:20]:
            issues.append(f"  '{rid}'")
        if len(bad_ids) > 20:
            issues.append(f"  ... and {len(bad_ids) - 20} more")
        issues.append("")

    # -- Required fields blank -------------------------------------------
    blank_issues: list[str] = []
    for row in rows:
        rid = row.get("resource_source_id", "(no id)")
        for field in REQUIRED_FIELDS:
            if not row.get(field, "").strip():
                blank_issues.append(f"  {rid}: '{field}' is blank")
    if blank_issues:
        issues.append("REQUIRED FIELDS WITH BLANK VALUES:")
        issues.extend(blank_issues[:50])
        if len(blank_issues) > 50:
            issues.append(f"  ... and {len(blank_issues) - 50} more")
        issues.append("")

    # -- Non-ASCII character detection -----------------------------------
    non_ascii_issues: list[str] = []
    for row in rows:
        rid = row.get("resource_source_id", "(no id)")
        for col, val in row.items():
            if not val:
                continue
            for i, ch in enumerate(val):
                if ord(ch) > 127:
                    suggestion = NON_ASCII_SUGGESTIONS.get(
                        ch,
                        f"(U+{ord(ch):04X} - replace with ASCII equivalent)"
                    )
                    context_start = max(0, i - 15)
                    context_end = min(len(val), i + 15)
                    context = val[context_start:context_end].replace("\t", " ")
                    non_ascii_issues.append(
                        f"  {rid} | {col} | char '{ch}' (U+{ord(ch):04X})\n"
                        f"    Context: ...{context}...\n"
                        f"    Suggestion: {suggestion}"
                    )
                    break  # one per cell is enough

    if non_ascii_issues:
        issues.append("NON-ASCII CHARACTERS DETECTED:")
        issues.append("  These may display as \ufffd (diamond question marks) or garbled text.")
        issues.append("  Fix them in the intermediate TSV before the next run.")
        issues.append("")
        issues.extend(non_ascii_issues[:40])
        if len(non_ascii_issues) > 40:
            issues.append(f"  ... and {len(non_ascii_issues) - 40} more cells with non-ASCII")
        issues.append("")

    return issues


def check_duplicate_uuids(rows: list[dict]) -> list[str]:
    """Check for duplicate UUIDs. Returns list of duplicate IDs if any.

    This is a critical error that should halt processing.
    """
    uuid_to_ids: dict[str, list[str]] = {}
    for row in rows:
        u = row.get("resource_uuid", "")
        rid = row.get("resource_source_id", "?")
        if u:
            uuid_to_ids.setdefault(u, []).append(rid)

    dupes = {u: ids for u, ids in uuid_to_ids.items() if len(ids) > 1}
    if dupes:
        messages = []
        for u, ids in dupes.items():
            messages.append(f"  UUID {u} -> source_ids: {', '.join(ids)}")
        return messages
    return []


def field_completion_summary(rows: list[dict]) -> list[str]:
    """Summarize field completion: populated vs blank for each column."""
    total = len(rows)
    lines: list[str] = []
    lines.append("FIELD COMPLETION")
    lines.append(f"  {'Field':<35s}  {'Populated':>9s}  {'Blank':>5s}")
    lines.append(f"  {'-' * 35}  {'-' * 9}  {'-' * 5}")
    for col in EXPECTED_COLUMNS:
        populated = sum(1 for r in rows if r.get(col, "").strip())
        blank = total - populated
        lines.append(f"  {col:<35s}  {populated:>9d}  {blank:>5d}")
    lines.append("")
    return lines


def facet_value_counts(rows: list[dict]) -> list[str]:
    """Count resources per value for faceted-filter fields.

    Parses semicolon-separated lists in resource_tool_type and
    resource_research_area, counting each row once per unique value.
    """
    lines: list[str] = []
    lines.append("EXPECTED FACET FILTER COUNTS")

    for field in ("resource_tool_type", "resource_research_area"):
        counts: Counter = Counter()
        for row in rows:
            val = row.get(field, "")
            for v in val.split(";"):
                v = v.strip()
                if v:
                    counts[v] += 1
        lines.append(f"\n  {field}  ({len(counts)} unique values):")
        for v, c in counts.most_common():
            lines.append(f"      {c:>4d}  {v}")

    lines.append("")
    return lines


def write_output(rows: list[dict], output_path: str) -> None:
    """Write processed rows to TSV with enforced column order."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=EXPECTED_COLUMNS, delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"[OK] Wrote {len(rows)} rows -> {output_path}")


def write_issues(issues: list[str], report_path: str, input_path: str,
                 rows: list[dict]) -> None:
    """Write the issue report to a text file."""
    os.makedirs(os.path.dirname(report_path), exist_ok=True)

    lines: list[str] = []
    lines.append("=" * 70)
    lines.append("  Resources Packaging Issue Report")
    lines.append(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"  Input file: {input_path}")
    lines.append(f"  Rows processed: {len(rows)}")
    lines.append("=" * 70)
    lines.append("")

    # -- Issues ----------------------------------------------------------
    if issues:
        lines.append("ISSUES")
        lines.append("-" * 70)
        lines.extend(issues)
    else:
        lines.append("No issues found. All checks passed.")

    lines.append("")
    lines.append("-" * 70)
    lines.append("")

    # -- Field completion ------------------------------------------------
    lines.extend(field_completion_summary(rows))

    # -- Facet filter counts ---------------------------------------------
    lines.extend(facet_value_counts(rows))

    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[OK] Issue report -> {report_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Orchestrate the resources packaging pipeline."""
    # 1. Find latest input file
    input_path = find_latest_input()
    if input_path is None:
        print(f"[FATAL] No resources_YYYY-MM-DD.tsv files found in {INPUT_DIR}")
        sys.exit(1)

    file_date = extract_date_from_path(input_path)
    print(f"[OK] Using input: {input_path} (date: {file_date})")

    # 2. Read and parse
    rows, extra_columns = read_resources(input_path)
    print(f"[OK] Read {len(rows)} rows")

    # 3. Fix type column
    type_issues = fix_type_column(rows)

    # 4. Generate UUIDs (overwrites any existing values)
    generate_uuids(rows)

    # 5. Check for duplicate UUIDs (critical)
    uuid_dupes = check_duplicate_uuids(rows)
    if uuid_dupes:
        print("[FATAL] Duplicate UUIDs detected (implies duplicate source IDs):")
        for msg in uuid_dupes:
            print(msg)
        sys.exit(1)

    # 6. Validate
    issues = validate(rows, extra_columns)
    if type_issues:
        issues.insert(0, "TYPE COLUMN AUTO-CORRECTIONS:")
        issues[1:1] = type_issues + [""]

    # 7. Write output
    output_path = os.path.join(OUTPUT_DIR, f"resources_{file_date}_processed.tsv")
    write_output(rows, output_path)

    # 8. Write issue report
    report_path = os.path.join(REPORT_DIR, f"resources_{file_date}_issue_report.txt")
    write_issues(issues, report_path, input_path, rows)

    # 9. Summary
    if issues:
        print(f"[WARN] {len(issues)} issue line(s) logged. Review: {report_path}")
    else:
        print("[OK] No issues found.")


if __name__ == "__main__":
    main()
