"""
test_package_resources.py
2026-05-14 ZD

Pytest test suite for the `package_resources.py` module.
"""

import csv
import os
import sys
import uuid

import pytest

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.package_resources import (
    EXPECTED_COLUMNS,
    ID_PATTERN,
    NAMESPACE,
    NON_ASCII_SUGGESTIONS,
    REQUIRED_FIELDS,
    check_duplicate_uuids,
    extract_date_from_path,
    facet_value_counts,
    field_completion_summary,
    find_latest_input,
    fix_type_column,
    generate_uuids,
    read_resources,
    sort_semicolon_lists,
    validate,
    write_issues,
    write_output,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_row(**overrides) -> dict:
    """Return a minimal valid row with optional field overrides."""
    defaults = {
        "type": "resource",
        "resource_uuid": "",
        "resource_source_id": "r4r_0001",
        "resource_title": "Test Resource",
        "resource_short_description": "A short description",
        "resource_source_url": "https://example.com",
        "resource_tool_type": "Analysis Tools",
        "resource_tool_subtype": "",
        "resource_research_area": "Cancer Biology",
        "resource_research_type": "",
        "resource_access": "Open",
        "resource_doc": "",
        "resource_poc_email": "",
        "resource_poc_name": "",
        "resource_full_description": "Full description here.",
    }
    defaults.update(overrides)
    return defaults


def _write_tsv(path: str, rows: list[dict],
               columns: list[str] | None = None) -> None:
    """Write rows to a TSV file at the given path."""
    cols = columns or EXPECTED_COLUMNS
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ===========================================================================
# find_latest_input
# ===========================================================================

class TestFindLatestInput:

    def test_returns_latest_date(self, tmp_path):
        """Picks the file with the most recent date in the filename."""
        for name in ("resources_2025-01-01.tsv", "resources_2026-03-15.tsv",
                      "resources_2026-01-10.tsv"):
            (tmp_path / name).write_text("header\n")

        result = find_latest_input(str(tmp_path))
        assert result is not None
        assert "2026-03-15" in os.path.basename(result)

    def test_ignores_non_matching_files(self, tmp_path):
        """Files that don't match resources_YYYY-MM-DD.tsv are skipped."""
        (tmp_path / "other_data.tsv").write_text("x\n")
        (tmp_path / "resources_latest.tsv").write_text("x\n")
        assert find_latest_input(str(tmp_path)) is None

    def test_empty_directory(self, tmp_path):
        assert find_latest_input(str(tmp_path)) is None


# ===========================================================================
# extract_date_from_path
# ===========================================================================

class TestExtractDateFromPath:

    def test_extracts_date(self):
        assert extract_date_from_path("data/resources_2026-05-14.tsv") == "2026-05-14"

    def test_raises_on_no_date(self):
        with pytest.raises(ValueError, match="Cannot extract date"):
            extract_date_from_path("data/resources_latest.tsv")


# ===========================================================================
# read_resources
# ===========================================================================

class TestReadResources:

    def test_reads_valid_tsv(self, tmp_path):
        """Reads rows, strips whitespace, and reports no extra columns."""
        row = _make_row(resource_title="  Padded Title  ")
        path = str(tmp_path / "resources_2026-01-01.tsv")
        _write_tsv(path, [row])

        rows, extra = read_resources(path)
        assert len(rows) == 1
        assert rows[0]["resource_title"] == "Padded Title"
        assert extra == []

    def test_detects_extra_columns(self, tmp_path):
        """Extra columns beyond EXPECTED_COLUMNS are reported."""
        row = _make_row()
        row["bonus_column"] = "surprise"
        path = str(tmp_path / "test.tsv")
        _write_tsv(path, [row], columns=EXPECTED_COLUMNS + ["bonus_column"])

        rows, extra = read_resources(path)
        assert "bonus_column" in extra
        # Extra column values are excluded from the returned rows
        assert "bonus_column" not in rows[0]

    def test_exits_on_missing_column(self, tmp_path):
        """Missing expected columns trigger SystemExit."""
        path = str(tmp_path / "bad.tsv")
        with open(path, "w") as fh:
            fh.write("resource_title\tresource_source_id\n")
            fh.write("Title\tid_001\n")

        with pytest.raises(SystemExit):
            read_resources(path)

    def test_reads_utf8_bom(self, tmp_path):
        """Handles UTF-8 BOM encoding (common from Excel)."""
        row = _make_row()
        path = str(tmp_path / "bom.tsv")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            writer = csv.DictWriter(fh, fieldnames=EXPECTED_COLUMNS,
                                    delimiter="\t")
            writer.writeheader()
            writer.writerow(row)

        rows, _ = read_resources(path)
        assert len(rows) == 1
        assert rows[0]["type"] == "resource"


# ===========================================================================
# generate_uuids
# ===========================================================================

class TestGenerateUuids:

    def test_deterministic(self):
        """Same source_id always produces the same UUID."""
        rows = [_make_row(resource_source_id="r4r_0042")]
        generate_uuids(rows)
        expected = str(uuid.uuid5(NAMESPACE, "r4r_0042"))
        assert rows[0]["resource_uuid"] == expected

    def test_different_ids_produce_different_uuids(self):
        rows = [_make_row(resource_source_id="r4r_0001"),
                _make_row(resource_source_id="r4r_0002")]
        generate_uuids(rows)
        assert rows[0]["resource_uuid"] != rows[1]["resource_uuid"]

    def test_empty_source_id_gets_blank_uuid(self):
        rows = [_make_row(resource_source_id="")]
        generate_uuids(rows)
        assert rows[0]["resource_uuid"] == ""

    def test_overwrites_existing_uuid(self):
        """generate_uuids replaces whatever was there before."""
        rows = [_make_row(resource_source_id="r4r_0001",
                          resource_uuid="old-garbage-uuid")]
        generate_uuids(rows)
        assert rows[0]["resource_uuid"] != "old-garbage-uuid"
        assert rows[0]["resource_uuid"] == str(uuid.uuid5(NAMESPACE, "r4r_0001"))


# ===========================================================================
# fix_type_column
# ===========================================================================

class TestFixTypeColumn:

    def test_correct_value_unchanged(self):
        rows = [_make_row(type="resource")]
        issues = fix_type_column(rows)
        assert rows[0]["type"] == "resource"
        assert issues == []

    def test_blank_gets_filled(self):
        rows = [_make_row(type="")]
        issues = fix_type_column(rows)
        assert rows[0]["type"] == "resource"
        assert issues == []  # blank fill is silent

    def test_wrong_value_corrected_and_reported(self):
        rows = [_make_row(type="dataset", resource_source_id="r4r_0099")]
        issues = fix_type_column(rows)
        assert rows[0]["type"] == "resource"
        assert len(issues) == 1
        assert "r4r_0099" in issues[0]
        assert "dataset" in issues[0]


# ===========================================================================
# sort_semicolon_lists
# ===========================================================================

class TestSortSemicolonLists:

    def test_already_sorted_no_issues(self):
        rows = [_make_row(resource_tool_type="Analysis Tools;Datasets and Databases")]
        issues = sort_semicolon_lists(rows)
        assert issues == []
        assert rows[0]["resource_tool_type"] == "Analysis Tools;Datasets and Databases"

    def test_unsorted_gets_sorted_and_reported(self):
        rows = [_make_row(
            resource_source_id="r4r_0010",
            resource_research_area="Cancer Treatment;Cancer Biology",
        )]
        issues = sort_semicolon_lists(rows)
        assert rows[0]["resource_research_area"] == "Cancer Biology;Cancer Treatment"
        assert len(issues) == 1
        assert "r4r_0010" in issues[0]
        assert "resource_research_area" in issues[0]

    def test_single_value_not_touched(self):
        """Fields without semicolons are left alone."""
        rows = [_make_row(resource_tool_type="Analysis Tools")]
        issues = sort_semicolon_lists(rows)
        assert issues == []

    def test_case_insensitive_sort(self):
        rows = [_make_row(resource_doc="Zebra;alpha")]
        issues = sort_semicolon_lists(rows)
        assert rows[0]["resource_doc"] == "alpha;Zebra"
        assert len(issues) == 1

    def test_multiple_fields_sorted_independently(self):
        rows = [_make_row(
            resource_tool_type="Lab Tools;Analysis Tools",
            resource_research_area="Cancer Treatment;Cancer Biology",
        )]
        issues = sort_semicolon_lists(rows)
        assert rows[0]["resource_tool_type"] == "Analysis Tools;Lab Tools"
        assert rows[0]["resource_research_area"] == "Cancer Biology;Cancer Treatment"
        assert len(issues) == 2

    def test_empty_field_ignored(self):
        rows = [_make_row(resource_tool_subtype="")]
        issues = sort_semicolon_lists(rows)
        assert issues == []

    def test_duplicates_removed_and_reported(self):
        rows = [_make_row(
            resource_source_id="r4r_0050",
            resource_research_area="Cancer Biology;Cancer Treatment;Cancer Biology",
        )]
        issues = sort_semicolon_lists(rows)
        assert rows[0]["resource_research_area"] == "Cancer Biology;Cancer Treatment"
        assert len(issues) == 1
        assert "r4r_0050" in issues[0]
        assert "duplicate" in issues[0].lower()

    def test_duplicates_removed_and_result_sorted(self):
        rows = [_make_row(
            resource_tool_type="Zebra;Alpha;Zebra;Alpha",
        )]
        issues = sort_semicolon_lists(rows)
        assert rows[0]["resource_tool_type"] == "Alpha;Zebra"
        assert len(issues) == 1
        assert "duplicate" in issues[0].lower()

    def test_dedup_only_no_sort_issue(self):
        """When duplicates are removed but remaining values are already sorted,
        only the dedup issue is reported (not a separate sort issue)."""
        rows = [_make_row(
            resource_research_area="Alpha;Alpha;Beta",
        )]
        issues = sort_semicolon_lists(rows)
        assert rows[0]["resource_research_area"] == "Alpha;Beta"
        assert len(issues) == 1
        assert "duplicate" in issues[0].lower()


# ===========================================================================
# validate
# ===========================================================================

class TestValidate:

    def test_clean_rows_no_issues(self):
        rows = [_make_row()]
        assert validate(rows, []) == []

    def test_extra_columns_reported(self):
        issues = validate([_make_row()], ["bonus_col"])
        assert any("EXTRA COLUMNS" in line for line in issues)
        assert any("bonus_col" in line for line in issues)

    def test_duplicate_source_ids_reported(self):
        rows = [_make_row(resource_source_id="r4r_0001"),
                _make_row(resource_source_id="r4r_0001")]
        issues = validate(rows, [])
        assert any("DUPLICATE" in line for line in issues)
        assert any("r4r_0001" in line for line in issues)

    def test_malformed_ids_reported(self):
        rows = [_make_row(resource_source_id="BAD-ID")]
        issues = validate(rows, [])
        assert any("FORMAT" in line for line in issues)
        assert any("BAD-ID" in line for line in issues)

    def test_valid_id_patterns_not_flagged(self):
        """IDs matching prefix_digits format are fine."""
        rows = [_make_row(resource_source_id="r4r_0001")]
        issues = validate(rows, [])
        assert not any("FORMAT" in line for line in issues)

    def test_blank_required_fields_reported(self):
        rows = [_make_row(resource_title="", resource_source_id="r4r_0010")]
        issues = validate(rows, [])
        text = "\n".join(issues)
        assert "REQUIRED FIELDS" in text
        assert "resource_title" in text
        assert "r4r_0010" in text

    def test_non_ascii_detected(self):
        rows = [_make_row(resource_title="Smart\u2019s Quotes")]
        issues = validate(rows, [])
        text = "\n".join(issues)
        assert "NON-ASCII" in text
        assert "U+2019" in text

    def test_non_ascii_suggestion_included(self):
        rows = [_make_row(resource_title="En\u2013dash")]
        issues = validate(rows, [])
        text = "\n".join(issues)
        assert "hyphen" in text.lower()

    def test_multiple_issue_types_combined(self):
        """Validate returns issues from multiple categories in one list."""
        rows = [_make_row(resource_source_id="BAD", resource_title="")]
        issues = validate(rows, ["extra"])
        text = "\n".join(issues)
        assert "EXTRA COLUMNS" in text
        assert "FORMAT" in text
        assert "REQUIRED FIELDS" in text


# ===========================================================================
# check_duplicate_uuids
# ===========================================================================

class TestCheckDuplicateUuids:

    def test_no_duplicates(self):
        rows = [_make_row(resource_uuid="aaa", resource_source_id="id_1"),
                _make_row(resource_uuid="bbb", resource_source_id="id_2")]
        assert check_duplicate_uuids(rows) == []

    def test_duplicates_detected(self):
        rows = [_make_row(resource_uuid="same", resource_source_id="id_1"),
                _make_row(resource_uuid="same", resource_source_id="id_2")]
        result = check_duplicate_uuids(rows)
        assert len(result) == 1
        assert "id_1" in result[0]
        assert "id_2" in result[0]


# ===========================================================================
# field_completion_summary
# ===========================================================================

class TestFieldCompletionSummary:

    def test_counts_populated_and_blank(self):
        rows = [
            _make_row(resource_poc_email="a@b.com"),
            _make_row(resource_poc_email=""),
            _make_row(resource_poc_email="c@d.com"),
        ]
        lines = field_completion_summary(rows)
        text = "\n".join(lines)

        assert "FIELD COMPLETION" in text
        # resource_poc_email: 2 populated, 1 blank
        # Find the line for resource_poc_email
        email_line = [l for l in lines if "resource_poc_email" in l
                      and "Field" not in l][0]
        assert "2" in email_line
        assert "1" in email_line

    def test_all_populated(self):
        rows = [_make_row(), _make_row()]
        lines = field_completion_summary(rows)
        # resource_title should show 2 populated, 0 blank
        title_line = [l for l in lines if "resource_title" in l
                      and "Field" not in l][0]
        assert "2" in title_line
        assert "0" in title_line


# ===========================================================================
# facet_value_counts
# ===========================================================================

class TestFacetValueCounts:

    def test_counts_single_values(self):
        rows = [
            _make_row(resource_tool_type="Analysis Tools"),
            _make_row(resource_tool_type="Analysis Tools"),
            _make_row(resource_tool_type="Datasets and Databases"),
        ]
        lines = facet_value_counts(rows)
        text = "\n".join(lines)
        assert "2" in text and "Analysis Tools" in text
        assert "1" in text and "Datasets and Databases" in text

    def test_parses_semicolon_separated(self):
        """Semicolon-separated values are counted individually."""
        rows = [_make_row(resource_research_area="Cancer Biology; Cancer Treatment")]
        lines = facet_value_counts(rows)
        text = "\n".join(lines)
        assert "Cancer Biology" in text
        assert "Cancer Treatment" in text

    def test_shows_unique_count(self):
        rows = [
            _make_row(resource_tool_type="A"),
            _make_row(resource_tool_type="B"),
            _make_row(resource_tool_type="C"),
        ]
        lines = facet_value_counts(rows)
        text = "\n".join(lines)
        assert "3 unique values" in text

    def test_empty_values_not_counted(self):
        rows = [_make_row(resource_tool_type="")]
        lines = facet_value_counts(rows)
        text = "\n".join(lines)
        assert "0 unique values" in text


# ===========================================================================
# write_output
# ===========================================================================

class TestWriteOutput:

    def test_writes_tsv_with_correct_columns(self, tmp_path):
        rows = [_make_row(), _make_row(resource_source_id="r4r_0002")]
        path = str(tmp_path / "out" / "output.tsv")
        write_output(rows, path)

        with open(path, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            written = list(reader)

        assert len(written) == 2
        assert list(reader.fieldnames) == EXPECTED_COLUMNS

    def test_excludes_extra_columns(self, tmp_path):
        """Extra keys in the row dicts are not written to the TSV."""
        row = _make_row()
        row["junk"] = "should not appear"
        path = str(tmp_path / "output.tsv")
        write_output([row], path)

        with open(path, encoding="utf-8") as fh:
            header = fh.readline()
        assert "junk" not in header


# ===========================================================================
# write_issues
# ===========================================================================

class TestWriteIssues:

    def test_report_with_no_issues(self, tmp_path):
        path = str(tmp_path / "report.txt")
        rows = [_make_row()]
        write_issues([], path, "input.tsv", rows)

        content = open(path, encoding="utf-8").read()
        assert "No issues found" in content
        assert "Rows processed: 1" in content
        assert "FIELD COMPLETION" in content
        assert "EXPECTED FACET FILTER COUNTS" in content

    def test_report_with_issues(self, tmp_path):
        path = str(tmp_path / "report.txt")
        rows = [_make_row()]
        write_issues(["Problem A", "Problem B"], path, "input.tsv", rows)

        content = open(path, encoding="utf-8").read()
        assert "ISSUES" in content
        assert "Problem A" in content
        assert "Problem B" in content
        assert "No issues found" not in content

    def test_report_includes_input_path(self, tmp_path):
        path = str(tmp_path / "report.txt")
        write_issues([], path, "data/special_input.tsv", [_make_row()])

        content = open(path, encoding="utf-8").read()
        assert "data/special_input.tsv" in content

    def test_row_count_matches_rows(self, tmp_path):
        path = str(tmp_path / "report.txt")
        rows = [_make_row() for _ in range(5)]
        write_issues([], path, "x.tsv", rows)

        content = open(path, encoding="utf-8").read()
        assert "Rows processed: 5" in content


# ===========================================================================
# ID_PATTERN
# ===========================================================================

class TestIdPattern:

    @pytest.mark.parametrize("valid_id", [
        "r4r_0001", "r4r_0299", "cbiit_catalog_001", "abc_123_456",
    ])
    def test_valid_ids_match(self, valid_id):
        assert ID_PATTERN.match(valid_id) is not None

    @pytest.mark.parametrize("invalid_id", [
        "R4R_001",       # uppercase
        "0001",          # no prefix
        "r4r",           # no underscore segment
        "r4r-001",       # hyphen instead of underscore
        "_r4r_001",      # leading underscore
        "r4r_",          # trailing underscore with nothing after
    ])
    def test_invalid_ids_rejected(self, invalid_id):
        assert ID_PATTERN.match(invalid_id) is None


# ===========================================================================
# Integration: end-to-end round-trip
# ===========================================================================

class TestEndToEnd:

    def test_full_pipeline(self, tmp_path):
        """Write a TSV, read it back, process, and verify output."""
        # Setup
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        output_path = str(tmp_path / "output.tsv")
        report_path = str(tmp_path / "report.txt")

        rows_in = [
            _make_row(resource_source_id="r4r_0001", type=""),
            _make_row(resource_source_id="r4r_0002",
                      resource_tool_type="Datasets and Databases"),
        ]
        input_path = str(input_dir / "resources_2026-05-14.tsv")
        _write_tsv(input_path, rows_in)

        # Find
        found = find_latest_input(str(input_dir))
        assert found == input_path

        # Read
        rows, extra = read_resources(input_path)
        assert len(rows) == 2

        # Fix type
        fix_type_column(rows)
        assert all(r["type"] == "resource" for r in rows)

        # UUIDs
        generate_uuids(rows)
        assert all(r["resource_uuid"] for r in rows)
        assert rows[0]["resource_uuid"] != rows[1]["resource_uuid"]

        # No duplicate UUIDs
        assert check_duplicate_uuids(rows) == []

        # Validate
        issues = validate(rows, extra)
        assert issues == []

        # Write
        write_output(rows, output_path)
        assert os.path.exists(output_path)

        # Report
        write_issues(issues, report_path, input_path, rows)
        content = open(report_path, encoding="utf-8").read()
        assert "No issues found" in content
        assert "FIELD COMPLETION" in content
        assert "EXPECTED FACET FILTER COUNTS" in content
        assert "Rows processed: 2" in content
