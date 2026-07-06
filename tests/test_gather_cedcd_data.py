"""
test_gather_cedcd_data.py
2026-07-06 ZD

Pytest test suite for the gather_cedcd_data.py module.
"""

import os
import sys

import pandas as pd
import pytest
from unittest.mock import patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.gather_cedcd_data import (
    get_newest_cohort_versions,
    clean_newlines,
    gather_cedcd_data,
)


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def sample_cedcd_df():
    """Minimal CEDCD cohort DataFrame with duplicate titles at different
    version IDs, mimicking real CEDCD export behavior."""
    return pd.DataFrame({
        "dataset_id": [100, 200, 300, 400],
        "dataset_title": ["Cohort A", "Cohort A", "Cohort B", "Cohort C"],
        "description": ["Old desc A", "New desc A", "Desc B", "Desc C"],
        "principal_investigators": [", Smith, John", "Doe, Jane",
                                    "Lee, Ann", "Park, Bo"],
        "primary_disease": ["Breast Cancer", "Breast Cancer",
                            "Lung Cancer", "Multiple"],
        "cohort_acronym": ["CA", "CA", "CB", "CC"],
        "participant_count": [1000, 1200, 500, 800],
    })


# ============================================================
# get_newest_cohort_versions
# ============================================================

class TestGetNewestCohortVersions:
    """Tests for filtering to newest version of each cohort."""

    def test_keeps_newest_by_dataset_id(self, sample_cedcd_df):
        result = get_newest_cohort_versions(sample_cedcd_df)

        # Should keep 3 unique titles
        assert len(result) == 3

        # For "Cohort A", should keep dataset_id 200 (newer), not 100
        cohort_a = result[result["dataset_title"] == "Cohort A"]
        assert len(cohort_a) == 1
        assert cohort_a.iloc[0]["dataset_id"] == 200

    def test_no_duplicates_returns_all(self):
        df = pd.DataFrame({
            "dataset_id": [1, 2, 3],
            "dataset_title": ["A", "B", "C"],
        })
        result = get_newest_cohort_versions(df)
        assert len(result) == 3

    def test_all_same_title_keeps_one(self):
        df = pd.DataFrame({
            "dataset_id": [10, 20, 30],
            "dataset_title": ["Same", "Same", "Same"],
        })
        result = get_newest_cohort_versions(df)
        assert len(result) == 1
        assert result.iloc[0]["dataset_id"] == 30


# ============================================================
# clean_newlines
# ============================================================

class TestCleanNewlines:
    """Tests for removing newline characters from DataFrame values."""

    def test_removes_standard_newlines(self):
        df = pd.DataFrame({"col": ["hello\nworld", "foo\r\nbar"]},
                          dtype="object")
        result = clean_newlines(df)
        assert result.iloc[0]["col"] == "hello world"
        assert result.iloc[1]["col"] == "foo bar"

    def test_removes_unicode_line_separators(self):
        df = pd.DataFrame({"col": ["line\u2028break", "para\u2029break"]},
                          dtype="object")
        result = clean_newlines(df)
        assert "\u2028" not in result.iloc[0]["col"]
        assert "\u2029" not in result.iloc[1]["col"]

    def test_removes_escaped_newlines(self):
        df = pd.DataFrame({"col": ["escaped\\nnewline"]}, dtype="object")
        result = clean_newlines(df)
        assert "\\n" not in result.iloc[0]["col"]

    def test_excludes_specified_columns(self):
        df = pd.DataFrame({
            "keep_clean": ["no\nnewline"],
            "leave_alone": ["keep\nnewline"],
        }, dtype="object")
        result = clean_newlines(df, exclude_cols=["leave_alone"])
        assert result.iloc[0]["keep_clean"] == "no newline"
        # Excluded column should still have the newline converted to string
        # but the \n removal is skipped
        assert "keep" in result.iloc[0]["leave_alone"]

    def test_numeric_columns_unaffected(self):
        df = pd.DataFrame({"num": [1, 2, 3], "text": ["a\nb", "c", "d"]})
        result = clean_newlines(df)
        assert result["num"].tolist() == [1, 2, 3]


# ============================================================
# gather_cedcd_data (integration)
# ============================================================

class TestGatherCedcdData:
    """Integration test for the main gather_cedcd_data function."""

    def test_produces_expected_output(self, sample_cedcd_df, tmp_path):
        """Verify the main function processes CEDCD data and writes output."""
        input_csv = tmp_path / "cedcd_input.csv"
        output_csv = tmp_path / "cedcd_output.csv"

        # Write the sample data as input CSV
        sample_cedcd_df.to_csv(input_csv, index=False)

        with patch("modules.gather_cedcd_data.config") as mock_config:
            mock_config.CEDCD_INPUT_CSV = str(input_csv)
            mock_config.CEDCD_INTERMED_CSV = str(output_csv)

            gather_cedcd_data()

        # Output should exist
        assert output_csv.exists()

        result = pd.read_csv(output_csv)

        # Should have 3 rows (deduplicated from 4)
        assert len(result) == 3

        # Should have expected hard-coded values
        assert all(result["type"] == "cedcd_dataset")
        assert all(result["dataset_source_repo"] == "CEDCD")

        # Should have URLs built from dataset_id
        assert all(result["dataset_source_url"].str.startswith(
            "https://cedcd.nci.nih.gov/cohort?id="))

        # Should have renamed primary_disease column
        assert "related_diseases" in result.columns

        # Leading comma in PI should be stripped
        cohort_a = result[result["dataset_title"] == "Cohort A"]
        if len(cohort_a) > 0:
            pi_val = str(cohort_a.iloc[0].get("principal_investigators", ""))
            assert not pi_val.startswith(",")

        # cohort_acronym should be dropped
        assert "cohort_acronym" not in result.columns
