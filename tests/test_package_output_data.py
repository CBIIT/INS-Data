"""
test_package_output_data.py
2026-07-06 ZD

Pytest test suite for the package_output_data.py module.
Focuses on the standardization and validation functions that finalize
all data for INS loading.
"""

import os

import pandas as pd
import pytest

from modules.package_output_data import (
    add_type_column,
    reorder_columns,
    validate_first_columns,
    replace_defined_characters,
    normalize_encoding,
    remove_html_tags,
    remove_html_tags_from_df,
    clean_html_entities,
    format_datetime_columns,
    process_special_characters,
    remove_nan_strings,
    validate_listlike_columns,
    validate_and_clean_unique_nodes,
    enforce_int_values,
    remove_publications_before_projects,
    standardize_data,
)


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def minimal_column_config():
    """Minimal column config matching the structure used in config.py."""
    return {
        "program": {
            "keep_and_rename": {
                "program_id": "program_id",
                "program_name": "program_name",
            },
            "node_id": "program_id",
            "list_like_cols": None,
            "datetime_cols": None,
            "html_tag_cols": None,
            "int_cols": None,
            "exclude_special_char_processing": [],
        },
        "test_type": {
            "keep_and_rename": {
                "type": "type",
                "node_id": "node_id",
                "link_id": "link_id",
                "title": "title",
                "value": "value",
            },
            "node_id": "node_id",
            "link_id": "link_id",
            "list_like_cols": ["value"],
            "datetime_cols": None,
            "html_tag_cols": ["title"],
            "int_cols": None,
            "exclude_special_char_processing": [],
        },
    }


@pytest.fixture
def sample_publications_df():
    """Publications DataFrame for early-publication removal tests."""
    return pd.DataFrame({
        "pmid": [111, 222, 333, 444],
        "coreproject": ["P1", "P1", "P2", "P2"],
        "title": ["Pub A", "Pub B", "Pub C", "Pub D"],
        "publication_date": ["2015-06-01", "2020-03-15",
                             "2018-01-01", "2022-06-01"],
    })


@pytest.fixture
def sample_projects_df():
    """Projects DataFrame for early-publication removal tests."""
    return pd.DataFrame({
        "project_id": ["P1", "P2"],
        "project_start_date": ["2019-01-01", "2020-01-01"],
    })


# ============================================================
# add_type_column
# ============================================================

class TestAddTypeColumn:
    """Tests for adding and mapping the type column."""

    def test_dataset_types_map_to_dataset(self):
        df = pd.DataFrame({"col": [1, 2]})
        for dtype in ("dbgap_dataset", "geo_dataset", "sra_dataset",
                       "cedcd_dataset", "ctd2_dataset", "dceg_dataset",
                       "nccr_dataset"):
            result = add_type_column(df, dtype)
            assert all(result["type"] == "dataset"), (
                f"{dtype} should map to 'dataset'")

    def test_non_dataset_types_keep_original(self):
        df = pd.DataFrame({"col": [1, 2]})
        for dtype in ("program", "grant", "project", "publication"):
            result = add_type_column(df, dtype)
            assert all(result["type"] == dtype)

    def test_does_not_modify_original(self):
        df = pd.DataFrame({"col": [1]})
        add_type_column(df, "program")
        assert "type" not in df.columns


# ============================================================
# reorder_columns
# ============================================================

class TestReorderColumns:
    """Tests for column reordering and validation."""

    def test_reorders_and_drops_extra_columns(self, minimal_column_config):
        df = pd.DataFrame({
            "program_id": ["P1"],
            "program_name": ["Test"],
            "extra_col": ["drop me"],
        })
        result = reorder_columns(df, minimal_column_config, "program")
        assert list(result.columns) == ["program_id", "program_name"]
        assert "extra_col" not in result.columns

    def test_missing_column_raises(self, minimal_column_config):
        df = pd.DataFrame({"program_id": ["P1"]})
        with pytest.raises(ValueError, match="Fields missing"):
            reorder_columns(df, minimal_column_config, "program")

    def test_invalid_datatype_raises(self, minimal_column_config):
        df = pd.DataFrame({"col": [1]})
        with pytest.raises(ValueError, match="Invalid datatype"):
            reorder_columns(df, minimal_column_config, "nonexistent")


# ============================================================
# validate_first_columns
# ============================================================

class TestValidateFirstColumns:
    """Tests for column ordering validation."""

    def test_correct_order_passes(self, minimal_column_config):
        df = pd.DataFrame({
            "type": ["test"],
            "node_id": ["N1"],
            "link_id": ["L1"],
            "title": ["T"],
            "value": ["V"],
        })
        # Should not raise
        validate_first_columns(df, minimal_column_config, "test_type")

    def test_wrong_order_raises(self, minimal_column_config):
        df = pd.DataFrame({
            "node_id": ["N1"],
            "type": ["test"],
            "link_id": ["L1"],
            "title": ["T"],
            "value": ["V"],
        })
        with pytest.raises(ValueError, match="not in the expected order"):
            validate_first_columns(df, minimal_column_config, "test_type")


# ============================================================
# replace_defined_characters
# ============================================================

class TestReplaceDefinedCharacters:
    """Tests for character replacement mapping."""

    def test_smart_quotes_replaced(self):
        assert replace_defined_characters("\u201cquoted\u201d") == "'quoted'"

    def test_em_dash_replaced(self):
        assert replace_defined_characters("a\u2014b") == "a-b"

    def test_superscripts_replaced(self):
        assert replace_defined_characters("10\u00B2") == "10^2"

    def test_subscripts_replaced(self):
        assert replace_defined_characters("H\u2082O") == "H2O"

    def test_nan_returns_nan(self):
        assert pd.isna(replace_defined_characters(pd.NA))

    def test_non_string_returns_as_is(self):
        assert replace_defined_characters(42) == 42


# ============================================================
# normalize_encoding
# ============================================================

class TestNormalizeEncoding:
    """Tests for encoding normalization."""

    def test_accented_characters_stripped(self):
        assert normalize_encoding("café") == "cafe"
        assert normalize_encoding("naïve") == "naive"

    def test_ascii_unchanged(self):
        assert normalize_encoding("hello world") == "hello world"

    def test_nan_returns_nan(self):
        assert pd.isna(normalize_encoding(pd.NA))


# ============================================================
# remove_html_tags
# ============================================================

class TestRemoveHtmlTags:
    """Tests for HTML tag removal."""

    def test_removes_italic_tags(self):
        assert remove_html_tags("<i>italic text</i>") == "italic text"

    def test_removes_bold_tags(self):
        assert remove_html_tags("<b>bold</b> normal") == "bold normal"

    def test_nested_tags(self):
        assert remove_html_tags("<b><i>nested</i></b>") == "nested"

    def test_no_tags_unchanged(self):
        assert remove_html_tags("plain text") == "plain text"


# ============================================================
# remove_nan_strings
# ============================================================

class TestRemoveNanStrings:
    """Tests for replacing 'nan' strings with blanks."""

    def test_removes_nan_string(self):
        df = pd.DataFrame({"col": ["nan", "NaN", "NAN", "real value", "banana"]})
        result = remove_nan_strings(df)
        assert result.iloc[0]["col"] == ""
        assert result.iloc[1]["col"] == ""
        assert result.iloc[2]["col"] == ""
        assert result.iloc[3]["col"] == "real value"
        assert result.iloc[4]["col"] == "banana"

    def test_preserves_numeric_nan(self):
        df = pd.DataFrame({"col": [1.0, float("nan"), 3.0]})
        result = remove_nan_strings(df)
        # Numeric NaN should not be converted to empty string
        assert pd.isna(result.iloc[1]["col"])


# ============================================================
# validate_listlike_columns
# ============================================================

class TestValidateListlikeColumns:
    """Tests for semicolon-separated list validation."""

    def test_removes_spaces_after_semicolons(self, minimal_column_config):
        df = pd.DataFrame({
            "type": ["test"],
            "node_id": ["N1"],
            "link_id": ["L1"],
            "title": ["T"],
            "value": ["A; B; C"],
        })
        result = validate_listlike_columns(
            df, minimal_column_config, "test_type")
        assert result.iloc[0]["value"] == "A;B;C"


# ============================================================
# validate_and_clean_unique_nodes
# ============================================================

class TestValidateAndCleanUniqueNodes:
    """Tests for duplicate node detection and cleanup."""

    def test_removes_true_duplicates(self, minimal_column_config, tmp_path):
        df = pd.DataFrame({
            "type": ["test", "test"],
            "node_id": ["N1", "N1"],
            "link_id": ["L1", "L1"],
            "title": ["Same", "Same"],
            "value": ["V", "V"],
        })
        from unittest.mock import patch
        with patch("modules.package_output_data.config") as mock_config:
            mock_config.REMOVED_DUPLICATES = str(tmp_path / "dups_")
            result = validate_and_clean_unique_nodes(
                df, minimal_column_config, "test_type")

        assert len(result) == 1

    def test_keeps_different_link_ids(self, minimal_column_config, tmp_path):
        """Rows with same node_id but different link_ids are valid
        many-to-many associations, not duplicates."""
        df = pd.DataFrame({
            "type": ["test", "test"],
            "node_id": ["N1", "N1"],
            "link_id": ["L1", "L2"],
            "title": ["Same", "Same"],
            "value": ["V", "V"],
        })
        from unittest.mock import patch
        with patch("modules.package_output_data.config") as mock_config:
            mock_config.REMOVED_DUPLICATES = str(tmp_path / "dups_")
            result = validate_and_clean_unique_nodes(
                df, minimal_column_config, "test_type")

        assert len(result) == 2

    def test_missing_node_id_column_raises(self, minimal_column_config):
        df = pd.DataFrame({"wrong_col": ["N1"]})
        with pytest.raises(ValueError, match="not found"):
            validate_and_clean_unique_nodes(
                df, minimal_column_config, "test_type")


# ============================================================
# enforce_int_values
# ============================================================

class TestEnforceIntValues:
    """Tests for integer type enforcement."""

    def test_converts_float_to_int(self):
        config = {"test": {"int_cols": ["count"]}}
        df = pd.DataFrame({"count": [1.0, 2.0, 3.0]})
        result = enforce_int_values(df, config, "test")
        assert result["count"].dtype.name == "Int64"

    def test_preserves_nan(self):
        config = {"test": {"int_cols": ["count"]}}
        df = pd.DataFrame({"count": [1.0, pd.NA, 3.0]})
        result = enforce_int_values(df, config, "test")
        assert pd.isna(result.iloc[1]["count"])

    def test_no_int_cols_config_skips(self):
        config = {"test": {"int_cols": None}}
        df = pd.DataFrame({"count": [1.5, 2.5]})
        result = enforce_int_values(df, config, "test")
        # Should be unchanged
        assert result.iloc[0]["count"] == 1.5


# ============================================================
# remove_publications_before_projects
# ============================================================

class TestRemovePublicationsBeforeProjects:
    """Tests for filtering publications published before project start dates.
    This is a critical data quality step for the final output."""

    def test_removes_early_publications(
            self, sample_publications_df, sample_projects_df, tmp_path):
        """Publications more than N days before the project start should
        be removed."""
        from unittest.mock import patch
        with patch("modules.package_output_data.config") as mock_config:
            mock_config.REMOVED_EARLY_PUBLICATIONS = str(
                tmp_path / "removed_early.csv")
            result = remove_publications_before_projects(
                sample_publications_df, sample_projects_df,
                day_diff_allowed=365)

        # Pub A (2015-06-01) is >365 days before P1 start (2019-01-01)
        assert 111 not in result["pmid"].values
        # Pub C (2018-01-01) is >365 days before P2 start (2020-01-01)
        assert 333 not in result["pmid"].values
        # Pub B and D should be kept
        assert 222 in result["pmid"].values
        assert 444 in result["pmid"].values

    def test_keeps_all_when_none(
            self, sample_publications_df, sample_projects_df):
        """Setting day_diff_allowed=None should skip filtering entirely."""
        result = remove_publications_before_projects(
            sample_publications_df, sample_projects_df,
            day_diff_allowed=None)
        assert len(result) == len(sample_publications_df)

    def test_keeps_all_when_zero(
            self, sample_publications_df, sample_projects_df):
        """Setting day_diff_allowed=0 should also skip (falsy value)."""
        result = remove_publications_before_projects(
            sample_publications_df, sample_projects_df,
            day_diff_allowed=0)
        assert len(result) == len(sample_publications_df)

    def test_saves_removed_report(
            self, sample_publications_df, sample_projects_df, tmp_path):
        """Removed publications should be saved to a report CSV."""
        report_path = tmp_path / "removed_early.csv"
        from unittest.mock import patch
        with patch("modules.package_output_data.config") as mock_config:
            mock_config.REMOVED_EARLY_PUBLICATIONS = str(report_path)
            remove_publications_before_projects(
                sample_publications_df, sample_projects_df,
                day_diff_allowed=365)
        assert report_path.exists()
        removed = pd.read_csv(report_path)
        assert len(removed) > 0

    def test_column_structure_preserved(
            self, sample_publications_df, sample_projects_df, tmp_path):
        """Output should have the same columns as input."""
        from unittest.mock import patch
        with patch("modules.package_output_data.config") as mock_config:
            mock_config.REMOVED_EARLY_PUBLICATIONS = str(
                tmp_path / "removed.csv")
            result = remove_publications_before_projects(
                sample_publications_df, sample_projects_df,
                day_diff_allowed=365)
        assert list(result.columns) == list(
            sample_publications_df.columns)

    def test_invalid_day_diff_type_raises(
            self, sample_publications_df, sample_projects_df):
        with pytest.raises(ValueError, match="Expected int"):
            remove_publications_before_projects(
                sample_publications_df, sample_projects_df,
                day_diff_allowed="365")


# ============================================================
# standardize_data (integration)
# ============================================================

class TestStandardizeData:
    """Integration test for the full standardization pipeline."""

    def test_full_pipeline_produces_clean_output(
            self, minimal_column_config, tmp_path):
        """Run standardize_data through all steps and verify the output
        meets loading requirements."""
        df = pd.DataFrame({
            "type": ["test_type"],
            "node_id": ["N1"],
            "link_id": ["L1"],
            "title": ["<i>Smart \u201cquotes\u201d</i> and caf\u00e9"],
            "value": ["A; B; C"],
        })

        from unittest.mock import patch
        with patch("modules.package_output_data.config") as mock_config:
            mock_config.REMOVED_DUPLICATES = str(tmp_path / "dups_")
            result = standardize_data(
                df, minimal_column_config, "test_type")

        # Type should be set
        assert result.iloc[0]["type"] == "test_type"

        # HTML tags should be removed from title
        assert "<i>" not in result.iloc[0]["title"]
        assert "</i>" not in result.iloc[0]["title"]

        # Smart quotes should be replaced
        assert "\u201c" not in result.iloc[0]["title"]
        assert "\u201d" not in result.iloc[0]["title"]

        # Semicolons in list-like fields should have no trailing spaces
        assert "; " not in result.iloc[0]["value"]

        # Columns should be in config order
        assert list(result.columns) == [
            "type", "node_id", "link_id", "title", "value"]

        # First columns should be type, node_id, link_id
        assert result.columns[0] == "type"
        assert result.columns[1] == "node_id"
        assert result.columns[2] == "link_id"


# ============================================================
# remove_html_tags_from_df
# ============================================================

class TestRemoveHtmlTagsFromDf:
    """Tests for HTML tag removal across DataFrame columns."""

    def test_removes_tags_in_configured_columns(self, minimal_column_config):
        df = pd.DataFrame({
            "type": ["test"],
            "node_id": ["N1"],
            "link_id": ["L1"],
            "title": ["<b>Bold Title</b>"],
            "value": ["<i>should stay</i>"],
        })
        result = remove_html_tags_from_df(
            df, minimal_column_config, "test_type")
        # title is in html_tag_cols
        assert result.iloc[0]["title"] == "Bold Title"
        # value is NOT in html_tag_cols — should be unchanged
        assert "<i>" in result.iloc[0]["value"]

    def test_missing_html_column_raises(self, minimal_column_config):
        df = pd.DataFrame({"wrong_col": ["<b>text</b>"]})
        with pytest.raises(ValueError, match="Expected html-tagged column"):
            remove_html_tags_from_df(df, minimal_column_config, "test_type")

    def test_no_html_cols_configured_skips(self):
        config = {"test": {"html_tag_cols": None}}
        df = pd.DataFrame({"col": ["<b>kept</b>"]})
        result = remove_html_tags_from_df(df, config, "test")
        assert "<b>" in result.iloc[0]["col"]


# ============================================================
# clean_html_entities
# ============================================================

class TestCleanHtmlEntities:
    """Tests for HTML entity decoding."""

    def test_decodes_numeric_entities(self):
        df = pd.DataFrame({"col": ["Smith&#8217;s Study"]}, dtype="object")
        result = clean_html_entities(df, "program")
        assert "\u2019" in result.iloc[0]["col"] or "'" in result.iloc[0]["col"]

    def test_decodes_named_entities(self):
        df = pd.DataFrame({"col": ["A &amp; B"]}, dtype="object")
        result = clean_html_entities(df, "program")
        assert "A & B" in result.iloc[0]["col"]

    def test_excludes_dbgap_description(self):
        df = pd.DataFrame({
            "description": ["Keep &amp; encoded"],
            "title": ["Decode &amp; this"],
        }, dtype="object")
        result = clean_html_entities(df, "dbgap_dataset")
        assert "&amp;" in result.iloc[0]["description"]
        assert "& " in result.iloc[0]["title"]


# ============================================================
# format_datetime_columns
# ============================================================

class TestFormatDatetimeColumns:
    """Tests for datetime column formatting."""

    def test_formats_datetime_to_yyyy_mm_dd(self):
        config = {"test": {"datetime_cols": ["date_col"]}}
        df = pd.DataFrame({"date_col": ["2023-06-15T12:00:00Z"]})
        result = format_datetime_columns(df, config, "test")
        assert result.iloc[0]["date_col"] == "2023-06-15"

    def test_handles_iso_date_strings(self):
        config = {"test": {"datetime_cols": ["date_col"]}}
        df = pd.DataFrame({"date_col": ["2023-01-01", "2023-06-15"]})
        result = format_datetime_columns(df, config, "test")
        assert result.iloc[0]["date_col"] == "2023-01-01"
        assert result.iloc[1]["date_col"] == "2023-06-15"

    def test_missing_column_raises(self):
        config = {"test": {"datetime_cols": ["missing_col"]}}
        df = pd.DataFrame({"other": ["2023-01-01"]})
        with pytest.raises(ValueError, match="Expected datetime column"):
            format_datetime_columns(df, config, "test")

    def test_none_config_skips(self):
        config = {"test": {"datetime_cols": None}}
        df = pd.DataFrame({"date_col": ["2023-01-01T12:00:00Z"]})
        result = format_datetime_columns(df, config, "test")
        # Should be unchanged
        assert result.iloc[0]["date_col"] == "2023-01-01T12:00:00Z"


# ============================================================
# process_special_characters
# ============================================================

class TestProcessSpecialCharacters:
    """Tests for the combined character replacement + encoding pipeline."""

    def test_processes_all_non_excluded_columns(self):
        config = {"test": {"exclude_special_char_processing": ["keep_col"]}}
        df = pd.DataFrame({
            "clean_col": ["caf\u00e9 \u201cquotes\u201d"],
            "keep_col": ["caf\u00e9 \u201cquotes\u201d"],
        })
        result = process_special_characters(df, config, "test")
        # clean_col should have accent stripped and quotes replaced
        assert "\u00e9" not in result.iloc[0]["clean_col"]
        assert "\u201c" not in result.iloc[0]["clean_col"]
        # keep_col should be untouched
        assert "\u00e9" in result.iloc[0]["keep_col"]
