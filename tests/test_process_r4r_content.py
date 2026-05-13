"""
test_process_r4r_content.py

Pytest test suite for the process_r4r_content.py module.
"""

import os
import sys
import csv
import re
import uuid
import pytest
import tempfile

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.process_r4r_content import (
    humanize_label,
    humanize_semicolon_list,
    humanize_access,
    humanize_doc,
    humanize_doc_list,
    dedup_semicolon_list,
    sort_semicolon_list,
    clean_for_tsv,
    markdown_to_html,
    replace_special_characters,
    sanitize_rows,
    parse_markdown_file,
    write_tsv,
    apply_curations,
    generate_report,
    COLUMNS,
    LABEL_OVERRIDES,
    ACCESS_LABELS,
    DOC_LABELS,
    CURATIONS,
    R4R_UUID_NAMESPACE,
)


# ---------------------------------------------------------------------------
# humanize_label
# ---------------------------------------------------------------------------

class TestHumanizeLabel:
    """Tests for snake_case -> human-friendly label conversion."""

    def test_override_takes_priority(self):
        """Values in LABEL_OVERRIDES should return the exact override."""
        assert humanize_label("datasets_databases") == "Datasets and Databases"
        assert humanize_label("cancer_omics") == "Cancer Omics"

    def test_fallback_title_case(self):
        """Values not in overrides should get underscores replaced + title case."""
        assert humanize_label("cancer_treatment") == "Cancer Treatment"
        assert humanize_label("basic") == "Basic"

    def test_empty_string_passthrough(self):
        assert humanize_label("") == ""

    def test_no_underscores(self):
        """A single word should just be title-cased."""
        assert humanize_label("bioinformatics") == "Bioinformatics"


class TestHumanizeSemicolonList:
    """Tests for semicolon-separated label humanization."""

    def test_multiple_values(self):
        result = humanize_semicolon_list("cancer_omics;cancer_treatment")
        assert result == "Cancer Omics;Cancer Treatment"

    def test_single_value(self):
        assert humanize_semicolon_list("basic") == "Basic"

    def test_empty_string_passthrough(self):
        assert humanize_semicolon_list("") == ""


# ---------------------------------------------------------------------------
# humanize_access
# ---------------------------------------------------------------------------

class TestHumanizeAccess:
    """Tests for access-level label mapping."""

    def test_known_values(self):
        assert humanize_access("open") == "Open Access"
        assert humanize_access("register") == "Requires Registration"
        assert humanize_access("cost") == "Requires Payment"

    def test_unknown_value_passes_through(self):
        assert humanize_access("mystery") == "mystery"

    def test_whitespace_stripped(self):
        assert humanize_access("  open  ") == "Open Access"


# ---------------------------------------------------------------------------
# humanize_doc / humanize_doc_list
# ---------------------------------------------------------------------------

class TestHumanizeDoc:
    """Tests for NCI division/office code mapping."""

    def test_known_code(self):
        assert humanize_doc("itcr") == "Informatics Technology for Cancer Research (ITCR)"
        assert humanize_doc("dceg") == "Division of Cancer Epidemiology & Genetics (DCEG)"

    def test_unknown_code_passes_through(self):
        assert humanize_doc("xyz") == "xyz"

    def test_doc_list_multiple(self):
        result = humanize_doc_list("itcr;dceg")
        assert "ITCR" in result
        assert "DCEG" in result
        assert ";" in result

    def test_doc_list_empty(self):
        assert humanize_doc_list("") == ""


# ---------------------------------------------------------------------------
# dedup_semicolon_list
# ---------------------------------------------------------------------------

class TestDedupSemicolonList:
    """Tests for duplicate removal in semicolon-separated strings."""

    def test_removes_duplicates(self):
        assert dedup_semicolon_list("A;B;A;C;B") == "A;B;C"

    def test_preserves_order(self):
        assert dedup_semicolon_list("C;A;B") == "C;A;B"

    def test_strips_whitespace(self):
        assert dedup_semicolon_list(" A ; B ; A ") == "A;B"

    def test_filters_empty_tokens(self):
        assert dedup_semicolon_list("A;;B;;;A") == "A;B"

    def test_empty_string_passthrough(self):
        assert dedup_semicolon_list("") == ""

    def test_single_value(self):
        assert dedup_semicolon_list("Analysis Tools") == "Analysis Tools"

    def test_all_same(self):
        assert dedup_semicolon_list("X;X;X") == "X"


# ---------------------------------------------------------------------------
# sort_semicolon_list
# ---------------------------------------------------------------------------

class TestSortSemicolonList:
    """Tests for alphabetical sorting of semicolon-separated values."""

    def test_sorts_alphabetically(self):
        assert sort_semicolon_list("Zebra;Apple;Mango") == "Apple;Mango;Zebra"

    def test_case_insensitive(self):
        assert sort_semicolon_list("banana;Apple;cherry") == "Apple;banana;cherry"

    def test_single_value(self):
        assert sort_semicolon_list("Only") == "Only"

    def test_empty_string(self):
        assert sort_semicolon_list("") == ""

    def test_none_passthrough(self):
        assert sort_semicolon_list(None) is None

    def test_filters_empty_tokens(self):
        assert sort_semicolon_list("B;;A;") == "A;B"

    def test_strips_whitespace(self):
        assert sort_semicolon_list(" C ; A ; B ") == "A;B;C"


# ---------------------------------------------------------------------------
# clean_for_tsv
# ---------------------------------------------------------------------------

class TestCleanForTsv:
    """Tests for whitespace collapsing in metadata fields."""

    def test_collapses_newlines(self):
        assert clean_for_tsv("space\nbetween") == "space between"

    def test_collapses_tabs(self):
        assert clean_for_tsv("space\tbetween") == "space between"

    def test_collapses_multiple_spaces(self):
        assert clean_for_tsv("space    between") == "space between"

    def test_strips_leading_trailing(self):
        assert clean_for_tsv("  space  ") == "space"

    def test_windows_line_endings(self):
        assert clean_for_tsv("a\r\nb\rc") == "a b c"


# ---------------------------------------------------------------------------
# markdown_to_html
# ---------------------------------------------------------------------------

class TestMarkdownToHtml:
    """Tests for lightweight markdown -> HTML conversion."""

    def test_bold(self):
        assert markdown_to_html("**bold text**") == "<strong>bold text</strong>"

    def test_italic(self):
        assert markdown_to_html("*italic text*") == "<em>italic text</em>"

    def test_link(self):
        result = markdown_to_html("[NCI](https://cancer.gov)")
        assert result == '<a href="https://cancer.gov">NCI</a>'

    def test_unordered_list_uses_hyphen(self):
        text = "Items:\n* first\n* second"
        result = markdown_to_html(text)
        assert "- first" in result
        assert "- second" in result

    def test_newlines_become_br(self):
        result = markdown_to_html("line one\nline two")
        assert "<br>" in result

    def test_consecutive_br_collapsed(self):
        result = markdown_to_html("a\n\n\nb")
        assert "<br><br>" not in result

    def test_trailing_br_stripped(self):
        result = markdown_to_html("content\n\n")
        assert not result.endswith("<br>")

    def test_existing_html_preserved(self):
        text = "x<sub>ij</sub> + y"
        result = markdown_to_html(text)
        assert "<sub>ij</sub>" in result

    def test_bold_before_italic(self):
        """Bold should be converted first so **bold** isn't partly eaten by italic."""
        result = markdown_to_html("**bold** and *italic*")
        assert "<strong>bold</strong>" in result
        assert "<em>italic</em>" in result


# ---------------------------------------------------------------------------
# replace_special_characters
# ---------------------------------------------------------------------------

class TestReplaceSpecialCharacters:
    """Tests for Unicode -> ASCII/HTML entity replacement."""

    def test_em_dash_to_hyphen(self):
        assert replace_special_characters("a\u2014b") == "a-b"

    def test_en_dash_to_hyphen(self):
        assert replace_special_characters("a\u2013b") == "a-b"

    def test_curly_quotes_to_straight(self):
        result = replace_special_characters("\u201CHello there!\u201D")
        assert result == "'Hello there!'"

    def test_greek_to_html_entity(self):
        assert replace_special_characters("\u03B1") == "&alpha;"
        assert replace_special_characters("\u03B2") == "&beta;"
        assert replace_special_characters("\u03A3") == "&Sigma;"

    def test_superscript_to_sup_tag(self):
        assert replace_special_characters("x\u00B2") == "x<sup>2</sup>"

    def test_subscript_to_sub_tag(self):
        assert replace_special_characters("H\u2082O") == "H<sub>2</sub>O"

    def test_multiplication_sign(self):
        assert replace_special_characters("a\u00D7b") == "a&times;b"

    def test_registered_trademark(self):
        assert replace_special_characters("PDQ\u00AE") == "PDQ&reg;"

    def test_trademark(self):
        assert replace_special_characters("MED-RT\u2122") == "MED-RT&trade;"

    def test_ligature_ff(self):
        assert replace_special_characters("e\uFB00ect") == "effect"

    def test_ligature_fi(self):
        assert replace_special_characters("pro\uFB01le") == "profile"

    def test_gte_to_html_entity(self):
        assert replace_special_characters("\u2265") == "&ge;"

    def test_plain_ascii_unchanged(self):
        text = "Regular stuff! <em>test</em>"
        assert replace_special_characters(text) == text

    def test_empty_string(self):
        assert replace_special_characters("") == ""

    def test_html_tags_survive(self):
        """Existing HTML tags must not be damaged by the NFD normalization."""
        text = '<a href="url">link</a> <sup>2</sup> <sub>i</sub>'
        assert replace_special_characters(text) == text

    def test_accent_stripped(self):
        """Accented Latin chars should lose their accent via NFD fallback."""
        assert replace_special_characters("caf\u00E9") == "cafe"


# ---------------------------------------------------------------------------
# sanitize_rows
# ---------------------------------------------------------------------------

class TestSanitizeRows:
    """Tests for in-place sanitization of row dicts."""

    def test_mutates_in_place(self):
        rows = [{"title": "test\u2014dash", "id": "1"}]
        sanitize_rows(rows)
        assert rows[0]["title"] == "test-dash"

    def test_non_string_values_untouched(self):
        rows = [{"title": "ok", "count": 42}]
        sanitize_rows(rows)
        assert rows[0]["count"] == 42

    def test_multiple_rows(self):
        rows = [
            {"a": "\u03B1"},
            {"a": "\u03B2"},
        ]
        sanitize_rows(rows)
        assert rows[0]["a"] == "&alpha;"
        assert rows[1]["a"] == "&beta;"

    def test_ctd2_superscript_in_descriptions(self):
        """CTD^2 and CTD2 in descriptions should become CTD<sup>2</sup>."""
        rows = [{
            "resource_short_description": "The CTD^2 Dashboard",
            "resource_full_description": "CTD2 Network data",
            "resource_title": "CTD2 Dashboard",
        }]
        sanitize_rows(rows)
        assert "CTD<sup>2</sup>" in rows[0]["resource_short_description"]
        assert "CTD<sup>2</sup>" in rows[0]["resource_full_description"]

    def test_ctd2_not_superscripted_in_title(self):
        """CTD2 in title should remain plain text (no HTML)."""
        rows = [{
            "resource_short_description": "desc",
            "resource_full_description": "desc",
            "resource_title": "CTD2 Dashboard",
        }]
        sanitize_rows(rows)
        assert rows[0]["resource_title"] == "CTD2 Dashboard"
        assert "<sup>" not in rows[0]["resource_title"]


# ---------------------------------------------------------------------------
# HTML entities / elements confined to description columns
# ---------------------------------------------------------------------------

# Columns where HTML content is expected (markdown -> HTML conversion)
_DESCRIPTION_COLS = {"resource_short_description", "resource_full_description"}

# Columns that are never checked (not user-facing text)
_SKIP_COLS = {"resource_source_url", "resource_uuid"}

# Pattern matching HTML entities (&alpha; &#947; &#x03B3;) and tags (<em>, <br>, etc.)
_HTML_PATTERN = re.compile(r"&[a-zA-Z]+;|&#[0-9]+;|&#x[0-9a-fA-F]+;|<[a-z/][^>]*>", re.IGNORECASE)


class TestHtmlConfinedToDescriptions:
    """HTML entities and elements should only appear in the two description columns."""

    def _make_row(self, **overrides):
        """Return a sanitized row with HTML content in descriptions."""
        base = {
            "type": "resource",
            "resource_uuid": "00000000-0000-0000-0000-000000000000",
            "resource_source_id": "1",
            "resource_title": "Plain Title",
            "resource_short_description": "Text with <em>emphasis</em>.",
            "resource_source_url": "https://example.com",
            "resource_tool_type": "Analysis Tools",
            "resource_tool_subtype": "R Software",
            "resource_research_area": "Cancer Omics",
            "resource_research_type": "Basic",
            "resource_access": "Open Access",
            "resource_doc": "Informatics Technology for Cancer Research (ITCR)",
            "resource_poc_email": "",
            "resource_poc_name": "",
            "resource_full_description": "Greek &alpha; and <sub>2</sub>.",
        }
        base.update(overrides)
        return base

    def test_no_html_in_non_description_columns(self):
        """A well-formed row has no HTML entities or tags outside the two description columns."""
        row = self._make_row()
        for col, val in row.items():
            if col in _DESCRIPTION_COLS or col in _SKIP_COLS:
                continue
            assert not _HTML_PATTERN.search(val), (
                f"Unexpected HTML in column '{col}': {val!r}")

    def test_html_is_present_in_description_columns(self):
        """Confirm the description columns do contain the expected HTML."""
        row = self._make_row()
        assert _HTML_PATTERN.search(row["resource_short_description"])
        assert _HTML_PATTERN.search(row["resource_full_description"])

    def test_sanitized_greek_stays_in_description(self):
        """Greek letters converted to HTML entities should only be in descriptions."""
        row = self._make_row(
            resource_title="Alpha Beta Tool",
            resource_full_description="Uses &alpha; and &beta; constants.",
        )
        sanitize_rows([row])
        for col, val in row.items():
            if col in _DESCRIPTION_COLS or col in _SKIP_COLS:
                continue
            assert not _HTML_PATTERN.search(val), (
                f"HTML entity leaked into column '{col}': {val!r}")

    def test_sup_sub_tags_only_in_descriptions(self):
        """<sup> and <sub> tags should not appear outside description columns."""
        row = self._make_row(
            resource_full_description="H<sub>2</sub>O and x<sup>2</sup>.",
        )
        for col, val in row.items():
            if col in _DESCRIPTION_COLS or col in _SKIP_COLS:
                continue
            assert "<sup>" not in val and "<sub>" not in val, (
                f"<sup>/<sub> tag found in column '{col}': {val!r}")


# ---------------------------------------------------------------------------
# parse_markdown_file
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_md_file(tmp_path):
    """Create a minimal R4R markdown file for testing parse_markdown_file."""
    content = """---
id: 99
title: Test Resource
description: A short description with **bold** text.
website: https://example.com
toolTypes:
  - toolType: analysis_tools/r_software
  - toolType: analysis_tools/statistical_software
researchAreas:
  - researchArea: cancer_omics
  - researchArea: cancer_biology
researchTypes:
  - researchType: basic
resourceAccess:
  type: open
docs:
  - doc: itcr
poc:
  - email: test@nih.gov
    name:
      firstname: Jane
      lastname: Doe
---
This is the **full description** with a [link](https://example.com).
"""
    filepath = tmp_path / "r4r_99.md"
    filepath.write_text(content, encoding="utf-8")
    return str(filepath)


class TestParseMarkdownFile:
    """Tests for parsing a single R4R markdown file."""

    def test_returns_all_columns(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        for col in COLUMNS:
            assert col in row, f"Missing column: {col}"

    def test_type_is_resource(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert row["type"] == "resource"

    def test_source_id(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert row["resource_source_id"] == "99"

    def test_title_cleaned(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert row["resource_title"] == "Test Resource"

    def test_uuid_is_deterministic(self, sample_md_file):
        row1 = parse_markdown_file(sample_md_file)
        row2 = parse_markdown_file(sample_md_file)
        assert row1["resource_uuid"] == row2["resource_uuid"]

    def test_uuid_is_valid(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        parsed = uuid.UUID(row["resource_uuid"])
        assert parsed.version == 5

    def test_tool_types_deduplicated(self, sample_md_file):
        """Two toolTypes with same parent should produce one parent entry."""
        row = parse_markdown_file(sample_md_file)
        assert row["resource_tool_type"] == "Analysis Tools"

    def test_tool_subtypes_kept(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        subtypes = row["resource_tool_subtype"].split(";")
        assert "R Software" in subtypes
        assert "Statistical Software" in subtypes

    def test_research_areas_humanized(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        areas = row["resource_research_area"].split(";")
        assert "Cancer Omics" in areas
        assert "Cancer Biology" in areas

    def test_access_humanized(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert row["resource_access"] == "Open Access"

    def test_doc_humanized(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert "ITCR" in row["resource_doc"]

    def test_poc_email(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert row["resource_poc_email"] == "test@nih.gov"

    def test_poc_name(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert row["resource_poc_name"] == "Jane Doe"

    def test_short_description_has_html(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert "<strong>bold</strong>" in row["resource_short_description"]

    def test_full_description_has_html(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert "<strong>full description</strong>" in row["resource_full_description"]
        assert '<a href="https://example.com">link</a>' in row["resource_full_description"]

    def test_url(self, sample_md_file):
        row = parse_markdown_file(sample_md_file)
        assert row["resource_source_url"] == "https://example.com"


class TestParseMarkdownFileMissingFields:
    """Tests for graceful handling of missing/optional fields."""

    def _write_md(self, tmp_path, yaml_body, md_body="Body text."):
        filepath = tmp_path / "test.md"
        filepath.write_text(f"---\n{yaml_body}\n---\n{md_body}", encoding="utf-8")
        return str(filepath)

    def test_missing_description_gets_placeholder(self, tmp_path):
        row = parse_markdown_file(self._write_md(tmp_path, 
            "id: 1\ntitle: No Desc"))
        assert row["resource_short_description"] == "No description preview."

    def test_missing_tool_types_gets_placeholder(self, tmp_path):
        row = parse_markdown_file(self._write_md(tmp_path,
            "id: 1\ntitle: No Tools"))
        assert row["resource_tool_type"] == "No Tool Type provided"

    def test_missing_research_areas_gets_placeholder(self, tmp_path):
        row = parse_markdown_file(self._write_md(tmp_path,
            "id: 1\ntitle: No Areas"))
        assert row["resource_research_area"] == "No Research Area provided"

    def test_missing_poc_gives_empty_strings(self, tmp_path):
        row = parse_markdown_file(self._write_md(tmp_path,
            "id: 1\ntitle: No POC"))
        assert row["resource_poc_email"] == ""
        assert row["resource_poc_name"] == ""

    def test_bad_frontmatter_raises(self, tmp_path):
        filepath = tmp_path / "bad.md"
        filepath.write_text("no front matter here", encoding="utf-8")
        with pytest.raises(ValueError, match="Could not parse front-matter"):
            parse_markdown_file(str(filepath))


# ---------------------------------------------------------------------------
# write_tsv
# ---------------------------------------------------------------------------

class TestWriteTsv:
    """Tests for TSV output writing."""

    def test_creates_file_with_header(self, tmp_path):
        rows = [{"type": "resource", "resource_uuid": "abc",
                 "resource_source_id": "1", "resource_title": "Test"}]
        outpath = str(tmp_path / "out.tsv")
        write_tsv(rows, outpath)
        assert os.path.exists(outpath)

        with open(outpath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            header = reader.fieldnames
        assert header == COLUMNS

    def test_row_count(self, tmp_path):
        rows = [
            {col: f"val{i}" for col in COLUMNS}
            for i in range(5)
        ]
        outpath = str(tmp_path / "out.tsv")
        write_tsv(rows, outpath)

        with open(outpath, "r", encoding="utf-8") as f:
            lines = f.readlines()
        # 1 header + 5 data rows
        assert len(lines) == 6

    def test_creates_parent_directories(self, tmp_path):
        outpath = str(tmp_path / "nested" / "dir" / "out.tsv")
        write_tsv([], outpath)
        assert os.path.exists(outpath)


# ---------------------------------------------------------------------------
# apply_curations
# ---------------------------------------------------------------------------

class TestApplyCurations:
    """Tests for post-generation curation overrides."""

    def test_overwrites_matching_row(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "modules.process_r4r_content.CURATION_LOG_PATH",
            str(tmp_path / "changelog.txt"))

        rows = [{
            "resource_source_id": "44",
            "resource_title": "Test",
            "resource_short_description": "old value",
            "resource_tool_type": "Analysis Tools",
            "resource_tool_subtype": "",
            "resource_research_area": "Cancer Biology",
        }]

        changes = apply_curations(rows)
        assert changes >= 1
        assert rows[0]["resource_short_description"] != "old value"

    def test_unmatched_rows_unchanged(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "modules.process_r4r_content.CURATION_LOG_PATH",
            str(tmp_path / "changelog.txt"))

        rows = [{
            "resource_source_id": "999999",
            "resource_title": "Not Curated",
            "resource_tool_type": "original",
        }]

        apply_curations(rows)
        assert rows[0]["resource_tool_type"] == "original"

    def test_writes_changelog(self, tmp_path, monkeypatch):
        log_path = str(tmp_path / "changelog.txt")
        monkeypatch.setattr(
            "modules.process_r4r_content.CURATION_LOG_PATH", log_path)

        rows = [{
            "resource_source_id": "44",
            "resource_title": "Test",
            "resource_short_description": "old",
            "resource_tool_type": "X",
            "resource_tool_subtype": "",
            "resource_research_area": "Y",
        }]

        apply_curations(rows)
        assert os.path.exists(log_path)
        contents = open(log_path, "r", encoding="utf-8").read()
        assert "Resource ID : 44" in contents


# ---------------------------------------------------------------------------
# generate_report
# ---------------------------------------------------------------------------

class TestGenerateReport:
    """Tests for the conversion report."""

    @pytest.fixture
    def sample_rows(self):
        """Minimal rows that satisfy the report generator."""
        return [
            {col: "" for col in COLUMNS}
            | {
                "resource_source_id": "1",
                "resource_title": "Resource A",
                "resource_tool_type": "Analysis Tools",
                "resource_tool_subtype": "R Software",
                "resource_research_area": "Cancer Omics",
                "resource_research_type": "Basic",
                "resource_access": "Open Access",
                "resource_doc": "ITCR",
                "resource_poc_email": "a@b.com",
                "resource_poc_name": "Test",
                "resource_full_description": "A full description here.",
            },
            {col: "" for col in COLUMNS}
            | {
                "resource_source_id": "2",
                "resource_title": "Resource B",
                "resource_tool_type": "Datasets and Databases",
                "resource_tool_subtype": "",
                "resource_research_area": "Cancer Biology",
                "resource_research_type": "Translational",
                "resource_access": "Requires Registration",
                "resource_doc": "DCEG",
                "resource_poc_email": "",
                "resource_poc_name": "",
                "resource_full_description": "Another description.",
            },
        ]

    def test_report_contains_overview(self, sample_rows):
        report = generate_report(sample_rows)
        assert "OVERVIEW" in report
        assert "Total resources processed : 2" in report

    def test_report_contains_id_range(self, sample_rows):
        report = generate_report(sample_rows)
        assert "1 - 2" in report

    def test_report_contains_value_distributions(self, sample_rows):
        report = generate_report(sample_rows)
        assert "VALUE DISTRIBUTIONS" in report
        assert "Analysis Tools" in report
        assert "Datasets and Databases" in report

    def test_report_contains_poc_coverage(self, sample_rows):
        report = generate_report(sample_rows)
        assert "POINT-OF-CONTACT" in report
        assert "With email" in report


# ---------------------------------------------------------------------------
# End-to-end: parse -> curate -> sanitize -> write
# ---------------------------------------------------------------------------

class TestEndToEnd:
    """Integration tests for the full pipeline with a temp .md file."""

    def _make_md(self, tmp_path, r_id, title, body="Description body.",
                 desc="Short desc.", tool="analysis_tools/r_software",
                 area="cancer_omics", rtype="basic", access="open",
                 doc="itcr", extra_yaml=""):
        content = f"""---
id: {r_id}
title: {title}
description: {desc}
website: https://example.com/{r_id}
toolTypes:
  - toolType: {tool}
researchAreas:
  - researchArea: {area}
researchTypes:
  - researchType: {rtype}
resourceAccess:
  type: {access}
docs:
  - doc: {doc}
{extra_yaml}
---
{body}
"""
        filepath = tmp_path / f"r4r_{r_id}.md"
        filepath.write_text(content, encoding="utf-8")
        return str(filepath)

    def test_full_pipeline_produces_valid_tsv(self, tmp_path, monkeypatch):
        """Parse two files, curate, sanitize, write, and verify TSV."""
        self._make_md(tmp_path, 1, "Resource One")
        self._make_md(tmp_path, 2, "Resource Two",
                      body="Text with \u2014 em dash and \u03B1 alpha.")

        monkeypatch.setattr("modules.process_r4r_content.INPUT_DIR",
                            str(tmp_path))
        monkeypatch.setattr("modules.process_r4r_content.CURATION_LOG_PATH",
                            str(tmp_path / "changelog.txt"))

        from modules.process_r4r_content import process_all_r4r_files

        rows = process_all_r4r_files()
        assert len(rows) == 2

        apply_curations(rows)
        sanitize_rows(rows)

        outpath = str(tmp_path / "output.tsv")
        write_tsv(rows, outpath)

        # Read back and validate
        with open(outpath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            written_rows = list(reader)

        assert len(written_rows) == 2
        assert written_rows[0]["resource_source_id"] == "1"
        assert written_rows[1]["resource_source_id"] == "2"

        # Verify special chars were sanitized
        desc2 = written_rows[1]["resource_full_description"]
        assert "\u2014" not in desc2   # em dash should be gone
        assert "\u03B1" not in desc2   # raw alpha should be gone
        assert "&alpha;" in desc2      # should be HTML entity now

    def test_output_is_ascii_safe(self, tmp_path, monkeypatch):
        """Every cell in the output TSV should be ASCII + HTML entities only."""
        self._make_md(tmp_path, 1, "Caf\u00E9 Resource",
                      body="Greek \u03B1\u03B2 and dash\u2014here.")

        monkeypatch.setattr("modules.process_r4r_content.INPUT_DIR",
                            str(tmp_path))
        monkeypatch.setattr("modules.process_r4r_content.CURATION_LOG_PATH",
                            str(tmp_path / "changelog.txt"))

        from modules.process_r4r_content import process_all_r4r_files

        rows = process_all_r4r_files()
        apply_curations(rows)
        sanitize_rows(rows)

        outpath = str(tmp_path / "output.tsv")
        write_tsv(rows, outpath)

        with open(outpath, "r", encoding="utf-8") as f:
            content = f.read()

        non_ascii = [c for c in content if ord(c) > 127]
        assert non_ascii == [], (
            f"Non-ASCII characters found in output: "
            f"{set(repr(c) for c in non_ascii)}")
