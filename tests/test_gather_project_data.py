"""
test_gather_project_data.py
2026-07-06 ZD

Pytest test suite for the gather_project_data.py module.
"""

import os
import pandas as pd
import pytest
from unittest.mock import patch

from modules.gather_project_data import (
    validate_identical_values,
    detect_supplement,
    drop_extra_supplement_rows,
    get_newest_or_oldest_value,
    get_list_values,
    gather_project_data,
)


@pytest.fixture
def sample_grants_df():
    """Minimal grants DataFrame matching the structure expected by
    gather_project_data."""
    return pd.DataFrame([
        {
            "queried_project_id": "R01CA111111",
            "grant_id": "R01CA111111-01",
            "project_title": "Old Title",
            "abstract_text": "Old abstract",
            "project_start_date": "2018-04-01",
            "project_end_date": "2020-03-31",
            "opportunity_number": "RFA-CA-18-001",
            "api_source_search": "nofo_RFA-CA-18-001",
            "org_name": "Test University",
            "org_city": "Bethesda",
            "org_state": "MD",
            "org_country": "UNITED STATES",
            "award_notice_date": "2018-01-15",
            "fiscal_year": 2018,
            "program.program_id": "PROG1",
        },
        {
            "queried_project_id": "R01CA111111",
            "grant_id": "R01CA111111-02",
            "project_title": "Newer Title",
            "abstract_text": "Newer abstract",
            "project_start_date": "2018-04-01",
            "project_end_date": "2021-03-31",
            "opportunity_number": "RFA-CA-18-001",
            "api_source_search": "nofo_RFA-CA-18-001",
            "org_name": "Test University",
            "org_city": "Bethesda",
            "org_state": "MD",
            "org_country": "UNITED STATES",
            "award_notice_date": "2019-01-15",
            "fiscal_year": 2019,
            "program.program_id": "PROG1",
        },
        {
            "queried_project_id": "U01CA222222",
            "grant_id": "U01CA222222-01",
            "project_title": "Second Project",
            "abstract_text": "Second abstract",
            "project_start_date": "2020-07-01",
            "project_end_date": "2023-06-30",
            "opportunity_number": "PAR-20-100",
            "api_source_search": "nofo_PAR-20-100",
            "org_name": "Other University",
            "org_city": "Rockville",
            "org_state": "MD",
            "org_country": "UNITED STATES",
            "award_notice_date": "2020-06-01",
            "fiscal_year": 2020,
            "program.program_id": "PROG2",
        },
    ])


# --- detect_supplement ---

class TestDetectSupplement:
    """Tests for supplement detection from grant IDs."""

    def test_standard_grant_not_supplement(self):
        assert detect_supplement("R01CA123456-01") is False

    def test_supplement_with_s(self):
        assert detect_supplement("3U24CA055727-26S1") is True

    def test_supplement_with_a(self):
        assert detect_supplement("1R01CA123456-01A1") is True

    def test_no_hyphen_with_a_in_id(self):
        # No hyphen means the whole string is checked; "A" in "CA" triggers True
        assert detect_supplement("R01CA123456") is True

    def test_no_hyphen_without_a_or_s(self):
        assert detect_supplement("R01HL123456") is False


# --- drop_extra_supplement_rows ---

class TestDropExtraSupplementRows:
    """Tests for dropping supplement rows when non-supplements exist."""

    def test_keeps_non_supplements(self):
        df = pd.DataFrame({
            "queried_project_id": ["P1", "P1", "P1"],
            "grant_id": ["P1-01", "P1-02", "P1-03S1"],
            "value": ["a", "b", "c"],
        })
        result = drop_extra_supplement_rows(df)
        assert len(result) == 2
        assert "P1-03S1" not in result["grant_id"].values

    def test_keeps_supplements_when_only_supplements(self):
        df = pd.DataFrame({
            "queried_project_id": ["P2", "P2"],
            "grant_id": ["P2-01S1", "P2-01A1"],
            "value": ["a", "b"],
        })
        result = drop_extra_supplement_rows(df)
        assert len(result) == 2


# --- validate_identical_values ---

class TestValidateIdenticalValues:
    """Tests for mismatch detection in fields expected to be identical."""

    def test_no_mismatch_returns_empty(self, sample_grants_df):
        result = validate_identical_values(sample_grants_df, "org_country")
        assert result is None or result.empty

    def test_mismatch_detected(self, sample_grants_df):
        df = sample_grants_df.copy()
        df.loc[0, "org_name"] = "Different University"
        result = validate_identical_values(df, "org_name")
        assert result is not None and not result.empty

    def test_skips_excluded_fields(self, sample_grants_df):
        assert validate_identical_values(
            sample_grants_df, "api_source_search") is None
        assert validate_identical_values(
            sample_grants_df, "program.program_id") is None


# --- get_newest_or_oldest_value ---

class TestGetNewestOrOldestValue:
    """Tests for retrieving newest/oldest values per project."""

    def test_newest_returns_latest_grant_value(self, sample_grants_df):
        result = get_newest_or_oldest_value(
            sample_grants_df, "project_title", "newest")
        assert result.loc["R01CA111111"] == "Newer Title"

    def test_oldest_returns_earliest_grant_value(self, sample_grants_df):
        result = get_newest_or_oldest_value(
            sample_grants_df, "project_title", "oldest")
        assert result.loc["R01CA111111"] == "Old Title"

    def test_invalid_agg_type_raises(self, sample_grants_df):
        with pytest.raises(ValueError, match="Invalid agg_type"):
            get_newest_or_oldest_value(
                sample_grants_df, "project_title", "invalid")


# --- get_list_values ---

class TestGetListValues:
    """Tests for collecting unique values as semicolon-separated strings."""

    def test_single_unique_value(self, sample_grants_df):
        result = get_list_values(sample_grants_df, "opportunity_number")
        assert result.loc["R01CA111111"] == "RFA-CA-18-001"

    def test_multiple_unique_values(self, sample_grants_df):
        df = sample_grants_df.copy()
        df.loc[1, "opportunity_number"] = "RFA-CA-19-002"
        result = get_list_values(df, "opportunity_number")
        assert "RFA-CA-18-001" in result.loc["R01CA111111"]
        assert "RFA-CA-19-002" in result.loc["R01CA111111"]


# --- gather_project_data (integration) ---

class TestGatherProjectData:
    """Integration test for the main gather_project_data function."""

    def test_produces_expected_output(self, sample_grants_df, tmp_path):
        """Verifies the main function aggregates grants into projects
        and writes output files."""
        output_csv = tmp_path / "project.csv"
        report_csv = tmp_path / "mismatch_report.csv"

        with patch("modules.gather_project_data.config") as mock_config:
            mock_config.PROJECTS_INTERMED_PATH = str(output_csv)
            mock_config.MISMATCHED_PROJECT_VALUES_REPORT = str(report_csv)

            result = gather_project_data(sample_grants_df)

        # Should produce one row per project-program combination
        assert len(result) == 2
        assert set(result["project_id"]) == {"R01CA111111", "U01CA222222"}

        # Newest title should be used
        r01_row = result[result["project_id"] == "R01CA111111"].iloc[0]
        assert r01_row["project_title"] == "Newer Title"

        # Output file should exist
        assert output_csv.exists()
