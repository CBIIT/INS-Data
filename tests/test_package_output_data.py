"""
test_package_output_data.py
2026-07-06 ZD

Pytest test suite for the package_output_data.py module.
Focuses on the standardization and validation functions that finalize
all data for INS loading.
"""

import os
import re

import pandas as pd
import pytest

import config as ins_config
import modules.package_output_data as package_module
from modules.package_output_data import (
    DATASET_TYPES,
    DATASET_HTML_ENTITY_EXCLUDE_COLS,
    ASCII_ONLY_COLS,
    read_csv_with_encoding_fallback,
    warn_on_encoding_damage,
    write_special_character_scan,
    is_preserved_character,
    add_type_column,
    reorder_columns,
    backfill_curated_dbgap_accessions,
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
    package_dbgap_datasets,
)


# Numeric (with or without trailing semicolon), hex, or named HTML entities
HTML_ENTITY_PATTERN = re.compile(
    r'&(?:#\d+;?|#[xX][0-9a-fA-F]+;?|[a-zA-Z]+;)')

# Entity-laden text modeled on real dbGaP API values
ENTITY_TEXT = ("Mej&#237;a &amp; Lo&#239;c&#8217;s &quot;study&quot; "
               "&#947;&#948;T&#160;cells")
# Final packaged value outside ASCII_ONLY_COLS: accented and Greek letters are
# kept; curly/straight double quotes and nbsp are mapped by
# replace_defined_characters
ENTITY_TEXT_DECODED = "Mej\u00eda & Lo\u00efc's 'study' \u03b3\u03b4T cells"


def build_raw_dataset_df(datatype):
    """Build a one-row pre-packaging DataFrame for a dataset datatype using
    the real config, filling every text field with HTML entities."""
    datatype_config = ins_config.COLUMN_CONFIGS[datatype]
    int_cols = datatype_config.get('int_cols') or []
    row = {}
    for raw_col, out_col in datatype_config['keep_and_rename'].items():
        if out_col in int_cols:
            row[raw_col] = [5]
        elif out_col == datatype_config['node_id']:
            row[raw_col] = ['D1']
        else:
            row[raw_col] = [ENTITY_TEXT]
    return pd.DataFrame(row)


def assert_no_disallowed_entities(df):
    """Fail if any HTML entity appears outside the allowed dataset fields."""
    for col in df.columns:
        if col in DATASET_HTML_ENTITY_EXCLUDE_COLS:
            continue
        for value in df[col].dropna().astype(str):
            assert not HTML_ENTITY_PATTERN.search(value), (
                f"HTML entity found in '{col}': {value!r}")


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
                       "nccr_dataset", "dctd_dataset"):
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


class TestBackfillCuratedDbgapAccessions:
    def test_adds_full_accession_from_gathered_data(self, tmp_path):
        curated = pd.DataFrame({
            "dataset_source_id": ["phs000001", "phs000002"]})
        gathered = pd.DataFrame({
            "accession": ["phs000001", "phs000002"],
            "full_accession": ["phs000001.v2.p1", "phs000002.v3.p2"],
        })
        gathered_path = tmp_path / "dbgap_datasets.csv"
        gathered.to_csv(gathered_path, index=False)

        result = backfill_curated_dbgap_accessions(
            curated, str(gathered_path))

        assert list(result["dataset_source_accession"]) == [
            "phs000001.v2.p1", "phs000002.v3.p2"]

    def test_preserves_existing_full_accessions(self):
        curated = pd.DataFrame({
            "dataset_source_accession": ["phs000001.v1.p1"]})

        result = backfill_curated_dbgap_accessions(
            curated, "/path/does/not/exist.csv")

        assert result.equals(curated)


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

    def test_preserve_keeps_greek_and_symbols(self):
        text = "\u03b2 \u03bcl 10\u00b2 H\u2082O \u2264 \u2265"
        assert replace_defined_characters(
            text, preserve_characters=True) == text

    def test_preserve_still_maps_quotes_dashes_and_spaces(self):
        result = replace_defined_characters(
            "\u201cA\u201d\u2014B\u00a0C \"D\" \u00a9", preserve_characters=True)
        assert result == "'A'-B C 'D' (C)"


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

    def test_decoded_entity_accent_handling_by_mode(self):
        """&#237; decodes to \u00ed, which replace_defined_characters never
        maps. normalize_encoding strips it to 'i' only in ASCII mode
        (ASCII_ONLY_COLS); in preserving mode (all other fields) it is kept."""
        decoded = "Mej\u00eda"
        assert replace_defined_characters(decoded) == decoded
        assert normalize_encoding(decoded) == "Mejia"
        assert normalize_encoding(decoded, preserve_characters=True) == decoded

    def test_preserve_composes_combining_accents(self):
        assert normalize_encoding(
            "Mej\u0069\u0301a", preserve_characters=True) == "Mej\u00eda"

    def test_preserve_keeps_allowed_and_folds_others(self):
        text = "Zo\u00eb 37\u00b0C 5\u00b12 \u03b3\u03b4T \u0141\u00f3d\u017a \u4e2d"
        assert normalize_encoding(text, preserve_characters=True) == (
            "Zo\u00eb 37\u00b0C 5\u00b12 \u03b3\u03b4T \u0141\u00f3d\u017a ")

    @pytest.mark.parametrize("text", [
        "Nguy\u1ec5n",                 # Latin Extended Additional
        "\u1e24asan",                  # Latin Extended Additional
        "\u1f04\u03bb\u03c6\u03b1",    # Greek Extended
    ])
    def test_preserve_keeps_extended_latin_and_greek(self, text):
        assert normalize_encoding(text, preserve_characters=True) == text

    def test_preserve_composes_extended_latin(self):
        assert normalize_encoding(
            "Nguye\u0302\u0303n", preserve_characters=True) == "Nguy\u1ec5n"


# ============================================================
# is_preserved_character
# ============================================================

class TestIsPreservedCharacter:
    """Tests for the allowlist of non-ASCII characters kept in output."""

    @pytest.mark.parametrize("char", list(
        "\u00ed\u00e9\u00f1\u00fc\u00f8\u0142\u0151\u00c9"
        "\u03b1\u03b2\u03b3\u03bc\u0394"
        "\u1ec5\u1e24\u1f04"
        "\u00b0\u00b2\u00b3\u00d7\u2264\u2265\u00b1\u00b5\u2082"))
    def test_preserved(self, char):
        assert is_preserved_character(char)

    @pytest.mark.parametrize("char", list(
        "\u201c\u2019\u2014\u00a0\u00ad\u0301\u00a9\u4e2d\u0439\ufffd\u1fbd"))
    def test_not_preserved(self, char):
        assert not is_preserved_character(char)

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

    def test_entities_decoded_before_char_normalization_and_tag_removal(
            self, minimal_column_config, monkeypatch):
        """clean_html_entities must run before process_special_characters
        (so decoded quotes are mapped like other text) and before
        remove_html_tags_from_df (so encoded tags are stripped)."""
        calls = []
        for name in ["clean_html_entities",
                     "process_special_characters",
                     "remove_html_tags_from_df"]:
            original = getattr(package_module, name)

            def spy(*args, _name=name, _original=original, **kwargs):
                calls.append(_name)
                return _original(*args, **kwargs)

            monkeypatch.setattr(package_module, name, spy)

        df = pd.DataFrame({
            "type": ["test_type"],
            "node_id": ["N1"],
            "link_id": ["L1"],
            "title": ["Smith&#8217;s &lt;i&gt;caf&#233;&lt;/i&gt;"],
            "value": ["Mej&#237;a"],
        })
        result = package_module.standardize_data(
            df, minimal_column_config, "test_type")

        assert calls.index("clean_html_entities") < calls.index(
            "process_special_characters")
        assert calls.index("clean_html_entities") < calls.index(
            "remove_html_tags_from_df")
        # Decoded ’ is only mapped to ' if decoding runs first; é/í are kept
        assert result.iloc[0]["title"] == "Smith's caf\u00e9"
        assert result.iloc[0]["value"] == "Mej\u00eda"

    @pytest.mark.parametrize("datatype", sorted(DATASET_TYPES))
    def test_dataset_outputs_have_no_entities_outside_allowed_fields(
            self, datatype):
        """Run the real config through standardize_data for every dataset
        type and confirm entities only survive in allowed fields."""
        df = build_raw_dataset_df(datatype)
        result = standardize_data(df, ins_config.COLUMN_CONFIGS, datatype)

        assert_no_disallowed_entities(result)

        datatype_config = ins_config.COLUMN_CONFIGS[datatype]
        int_cols = datatype_config.get('int_cols') or []
        skip_cols = (set(DATASET_HTML_ENTITY_EXCLUDE_COLS)
                     | set(int_cols)
                     | set(datatype_config.get(
                         'exclude_special_char_processing', []))
                     | {'type', datatype_config['node_id']})
        for col in result.columns:
            if col not in skip_cols:
                assert result.iloc[0][col] == ENTITY_TEXT_DECODED, col

        for col in DATASET_HTML_ENTITY_EXCLUDE_COLS:
            if col in result.columns:
                assert "&#237;" in result.iloc[0][col]

    def test_nullable_integer_columns_serialize_without_decimals(self):
        config = {
            "dataset": {
                "keep_and_rename": {
                    "type": "type",
                    "dataset_uuid": "dataset_uuid",
                    "participant_count": "participant_count",
                    "enrollment_year": "enrollment_year",
                },
                "node_id": "dataset_uuid",
                "list_like_cols": None,
                "datetime_cols": None,
                "html_tag_cols": None,
                "int_cols": ["participant_count", "enrollment_year"],
                "exclude_special_char_processing": [],
            }
        }
        df = pd.DataFrame({
            "type": ["dataset", "dataset"],
            "dataset_uuid": ["D1", "D2"],
            "participant_count": [89656.0, float("nan")],
            "enrollment_year": [1993.0, float("nan")],
        })

        result = standardize_data(df, config, "dataset")
        serialized = result.to_csv(sep="\t", index=False)

        assert result["participant_count"].dtype.name == "Int64"
        assert result["enrollment_year"].dtype.name == "Int64"
        assert "89656.0" not in serialized
        assert "1993.0" not in serialized
        assert "89656\t1993" in serialized


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

    def test_decodes_default_pandas_string_dtype(self):
        """Regression: pandas 3 text columns are StringDtype, not object."""
        df = pd.DataFrame({"PI_name": ["Juan Mej&#237;a-Arangure"]})
        assert df["PI_name"].dtype != object
        result = clean_html_entities(df, "dbgap_dataset")
        # process_special_characters keeps í outside ASCII_ONLY_COLS
        assert result.iloc[0]["PI_name"] == "Juan Mej\u00eda-Arangure"

    def test_decodes_entities_without_semicolon(self):
        df = pd.DataFrame({"col": ["Mej&#237a"]})
        result = clean_html_entities(df, "geo_dataset")
        assert result.iloc[0]["col"] == "Mej\u00eda"

    def test_preserves_missing_values(self):
        df = pd.DataFrame({"col": ["Men&#39;s", None]})
        result = clean_html_entities(df, "dbgap_dataset")
        assert result.iloc[0]["col"] == "Men's"
        assert pd.isna(result.iloc[1]["col"])

    @pytest.mark.parametrize("datatype", sorted(DATASET_TYPES))
    def test_dataset_types_retain_entities_only_in_allowed_fields(
            self, datatype):
        df = pd.DataFrame({
            "description": ["Keep &#39; encoded"],
            "experimental_approaches": ["Keep &#39; encoded"],
            "dataset_title": ["Men&#39;s Study"],
        })
        result = clean_html_entities(df, datatype)
        assert result.iloc[0]["description"] == "Keep &#39; encoded"
        assert result.iloc[0]["experimental_approaches"] == "Keep &#39; encoded"
        assert result.iloc[0]["dataset_title"] == "Men's Study"

    def test_non_dataset_types_decode_description(self):
        df = pd.DataFrame({"description": ["Men&#39;s"]})
        result = clean_html_entities(df, "project")
        assert result.iloc[0]["description"] == "Men's"


# ============================================================
# package_dbgap_datasets (end-to-end TSV output)
# ============================================================

class TestPackageDbgapDatasetsEntities:
    """End-to-end check that HTML entities from the dbGaP API do not reach
    the packaged TSV outside of allowed fields."""

    @pytest.mark.parametrize("dbgap_curated", [False, True])
    def test_written_tsv_has_no_disallowed_entities(
            self, tmp_path, dbgap_curated):
        from unittest.mock import patch

        # Round-trip through CSV to mirror how intermediate data is loaded
        input_path = tmp_path / "dbgap_intermediate.csv"
        build_raw_dataset_df("dbgap_dataset").to_csv(input_path, index=False)
        df = read_csv_with_encoding_fallback(input_path)

        output_path = tmp_path / "dbgap_datasets.tsv"
        with patch("modules.package_output_data.config") as mock_config:
            mock_config.DBGAP_OUTPUT_PATH = str(output_path)
            mock_config.DBGAP_OUTPUT_CURATED_CLEANED = str(output_path)
            package_dbgap_datasets(
                df, ins_config.COLUMN_CONFIGS, dbgap_curated=dbgap_curated)

        written = pd.read_csv(output_path, sep="\t", dtype=str,
                              keep_default_na=False, encoding="utf-8")
        assert_no_disallowed_entities(written)
        assert written.iloc[0]["PI_name"] == ENTITY_TEXT_DECODED
        assert "&#237;" in written.iloc[0]["description"]
        # Loader TSVs stay plain UTF-8 without a BOM
        assert not output_path.read_bytes().startswith(b"\xef\xbb\xbf")


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
        # clean_col keeps the accent but has quotes replaced
        assert result.iloc[0]["clean_col"] == "caf\u00e9 'quotes'"
        # keep_col should be untouched
        assert result.iloc[0]["keep_col"] == "caf\u00e9 \u201cquotes\u201d"

    @pytest.mark.parametrize("ascii_col", ASCII_ONLY_COLS)
    def test_ascii_only_cols_keep_legacy_processing(self, ascii_col):
        config = {"test": {"exclude_special_char_processing": []}}
        text = "Mej\u00eda 37\u00b0C \u03b3 \u2264 10\u00b2"
        df = pd.DataFrame({ascii_col: [text], "PI_name": [text]})
        result = process_special_characters(df, config, "test")
        assert result.iloc[0][ascii_col] == "Mejia 37C g <= 10^2"
        assert result.iloc[0]["PI_name"] == text


# ============================================================
# read_csv_with_encoding_fallback
# ============================================================

class TestReadCsvWithEncodingFallback:
    """Tests for the utf-8/cp1252 fallback used on Excel-edited curated files."""

    def test_reads_valid_utf8_file(self, tmp_path):
        path = tmp_path / "utf8.tsv"
        path.write_text("col_a\tcol_b\nfoo\tbar\n", encoding="utf-8")
        result = read_csv_with_encoding_fallback(path, sep="\t")
        assert result.iloc[0]["col_a"] == "foo"

    def test_falls_back_to_cp1252_on_decode_error(self, tmp_path):
        """cp1252-only bytes like 0xa0 (nbsp) and 0x96 (en dash), as commonly
        produced by Excel, raise UnicodeDecodeError under utf-8."""
        path = tmp_path / "cp1252.tsv"
        content = b"col_a\tcol_b\nfoo\xa0bar\tvalue\x96dash\n"
        path.write_bytes(content)

        with pytest.raises(UnicodeDecodeError):
            pd.read_csv(path, sep="\t", encoding="utf-8")

        result = read_csv_with_encoding_fallback(path, sep="\t")
        assert result.iloc[0]["col_a"] == "foo\xa0bar"
        assert result.iloc[0]["col_b"] == "value\u2013dash"

    def test_passes_through_kwargs(self, tmp_path):
        path = tmp_path / "dtype.tsv"
        path.write_text("id\tvalue\n1\t007\n", encoding="utf-8")
        result = read_csv_with_encoding_fallback(
            path, sep="\t", dtype={"value": str})
        assert result.iloc[0]["value"] == "007"

    def test_strips_utf8_bom(self, tmp_path):
        path = tmp_path / "bom.tsv"
        path.write_text("type\tPI_name\ndataset\tMej\u00eda\n",
                        encoding="utf-8-sig")
        result = read_csv_with_encoding_fallback(path, sep="\t")
        assert list(result.columns) == ["type", "PI_name"]
        assert result.iloc[0]["PI_name"] == "Mej\u00eda"

    def test_reads_excel_unicode_text_utf16(self, tmp_path):
        """Excel 'Unicode Text' saves UTF-16 and keeps every character."""
        path = tmp_path / "unicode_text.txt"
        path.write_text("PI_name\tterm\nMej\u00eda\tROR\u03b3t \u2264 5\n",
                        encoding="utf-16")
        result = read_csv_with_encoding_fallback(path, sep="\t")
        assert result.iloc[0]["PI_name"] == "Mej\u00eda"
        assert result.iloc[0]["term"] == "ROR\u03b3t \u2264 5"

    def test_excel_ansi_save_round_trip(self, tmp_path, capsys):
        """Simulates Excel 'Text (Tab delimited)': cp1252 keeps accented
        letters but turns characters like \u03b3 into '?', which is flagged."""
        original = pd.DataFrame({"PI_name": ["Mej\u00eda"],
                                 "term": ["ROR\u03b3t at 37\u00b0C"]})
        text = original.to_csv(sep="\t", index=False)
        path = tmp_path / "excel_ansi.tsv"
        path.write_bytes(text.encode("cp1252", errors="replace"))

        result = read_csv_with_encoding_fallback(path, sep="\t")

        assert result.iloc[0]["PI_name"] == "Mej\u00eda"
        assert result.iloc[0]["term"] == "ROR?t at 37\u00b0C"
        output = capsys.readouterr().out
        assert "Retrying with cp1252" in output
        assert "Excel row 2, term" in output


# ============================================================
# warn_on_encoding_damage
# ============================================================

class TestWarnOnEncodingDamage:
    """Tests for detection of likely Excel encoding damage."""

    def test_flags_question_mark_substitution(self):
        df = pd.DataFrame({"desc": ["fine", "Th17 ROR?t+ cells"]})
        findings = warn_on_encoding_damage(df, "test.tsv")
        assert [(f["source_file"], f["excel_row"], f["column"])
                for f in findings] == [("test.tsv", 3, "desc")]

    def test_flags_mojibake_and_replacement_char(self):
        df = pd.DataFrame({"name": ["Mej\u00c3\u00ada", "It\u00e2\u20ac\u2122s",
                                    "Mej\ufffda"]})
        findings = warn_on_encoding_damage(df, "test.tsv")
        assert sorted(f["excel_row"] for f in findings) == [2, 3, 4]

    def test_ignores_urls_and_clean_text(self):
        df = pd.DataFrame({
            "url": ["https://example.org/study.cgi?study_id=phs000001"],
            "html": ["<a href='./study.cgi?study_id=phs0'>link</a>"],
            "text": ["Mej\u00eda \u00c9cole RORγt \u2264 5? Yes"],
        })
        assert warn_on_encoding_damage(df, "test.tsv") == []

    @pytest.mark.parametrize("value", [
        "<i>NF?B</i> signaling",
        "NF?B/RelA complex",
        "Ihhi?EC/?iEC cells",
        "<a href='https://example.org/a?b=1'>NF?B</a>",
        "See www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1 for NF?B",
    ])
    def test_flags_damage_next_to_slashes_html_and_urls(self, value):
        findings = warn_on_encoding_damage(
            pd.DataFrame({"title": [value]}), "test.tsv")
        assert len(findings) == 1
        assert "NF?B" in findings[0]["snippet"] or "i?E" in findings[0]["snippet"]

    def test_console_truncated_but_all_findings_returned(self, capsys):
        df = pd.DataFrame({"desc": [f"NF?B row {i}" for i in range(15)]})
        findings = warn_on_encoding_damage(df, "test.tsv", max_examples=10)
        assert len(findings) == 15
        assert capsys.readouterr().out.count("Excel row") == 10


# ============================================================
# write_special_character_scan
# ============================================================

class TestWriteSpecialCharacterScan:
    """Tests for the full special character CSV report."""

    def test_collects_all_findings_across_files(self, tmp_path):
        findings = []
        for name, rows in [("dbgap.tsv", 12), ("sra.tsv", 3)]:
            path = tmp_path / name
            pd.DataFrame({"title": [f"NF?B {i}" for i in range(rows)]}).to_csv(
                path, sep="\t", index=False)
            read_csv_with_encoding_fallback(path, sep="\t", findings=findings)

        report_path = tmp_path / "reports" / "special_character_scan.csv"
        write_special_character_scan(findings, str(report_path))

        report = pd.read_csv(report_path, encoding="utf-8-sig")
        assert list(report.columns) == [
            "source_file", "issue", "excel_row", "column", "snippet"]
        assert len(report) == 15
        assert report["source_file"].str.endswith("dbgap.tsv").sum() == 12
        assert report["source_file"].str.endswith("sra.tsv").sum() == 3

    def test_writes_headers_when_clean(self, tmp_path):
        report_path = tmp_path / "special_character_scan.csv"
        write_special_character_scan([], str(report_path))
        report = pd.read_csv(report_path, encoding="utf-8-sig")
        assert report.empty
        assert "source_file" in report.columns
