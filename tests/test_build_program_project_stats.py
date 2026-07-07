"""
test_build_program_project_stats.py
2026-07-07 ZD

Pytest test suite for the build_program_project_stats.py module.
"""

import os
import sys

import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.build_program_project_stats import (
    get_grant_stats_by_program,
    get_shared_projects_by_program_pair,
    build_program_project_stats,
)


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def sample_grants():
    """Sample grants DataFrame with multiple programs and shared projects."""
    return pd.DataFrame({
        "program.program_id": ["PROG_A", "PROG_A", "PROG_A",
                               "PROG_B", "PROG_B",
                               "PROG_C"],
        "api_source_search": ["nofo_RFA-CA-01-001", "nofo_RFA-CA-01-001",
                              "award_R01CA111111", "nofo_RFA-CA-02-002",
                              "nofo_RFA-CA-02-002", "nofo_RFA-CA-03-003"],
        "queried_project_id": ["R01CA111111", "R01CA222222", "R01CA111111",
                               "R01CA222222", "R01CA333333", "R01CA444444"],
        "grant_id": ["R01CA111111-01", "R01CA222222-01", "R01CA111111-02",
                     "R01CA222222-02", "R01CA333333-01", "R01CA444444-01"],
        "fiscal_year": [2015, 2018, 2016, 2019, 2020, 2021],
    })


@pytest.fixture
def single_program_grants():
    """Grants DataFrame with only one program (no shared projects possible)."""
    return pd.DataFrame({
        "program.program_id": ["PROG_ONLY", "PROG_ONLY", "PROG_ONLY"],
        "api_source_search": ["nofo_1", "nofo_1", "nofo_2"],
        "queried_project_id": ["PRJ_X", "PRJ_Y", "PRJ_X"],
        "grant_id": ["G1", "G2", "G3"],
        "fiscal_year": [2020, 2021, 2022],
    })


# ============================================================
# get_grant_stats_by_program
# ============================================================

class TestGetGrantStatsByProgram:
    """Tests for aggregating grant statistics by program."""

    def test_output_csv_is_written(self, sample_grants, tmp_path):
        """Verify that the report CSV is created."""
        output_file = tmp_path / "grantsStatsByProgram.csv"
        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_AGG_FUNCS_BY_COL = {
                "api_source_search": "nunique",
                "queried_project_id": "nunique",
                "grant_id": "nunique",
                "fiscal_year": "min",
            }
            mock_cfg.STAT_FISCALYEAR_COL = "fiscal_year"
            mock_cfg.STAT_GRANTS_BY_PROGRAM_FILENAME = str(output_file)

            get_grant_stats_by_program(sample_grants)

        assert output_file.exists()

    def test_aggregated_values_are_correct(self, sample_grants, tmp_path):
        """Verify aggregation counts and earliest fiscal year per program."""
        output_file = tmp_path / "grantsStatsByProgram.csv"
        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_AGG_FUNCS_BY_COL = {
                "api_source_search": "nunique",
                "queried_project_id": "nunique",
                "grant_id": "nunique",
                "fiscal_year": "min",
            }
            mock_cfg.STAT_FISCALYEAR_COL = "fiscal_year"
            mock_cfg.STAT_GRANTS_BY_PROGRAM_FILENAME = str(output_file)

            get_grant_stats_by_program(sample_grants)

        result = pd.read_csv(output_file)

        # PROG_A: 2 unique searches, 2 unique projects, 3 unique grants, min FY 2015
        prog_a = result[result["program.program_id"] == "PROG_A"].iloc[0]
        assert prog_a["api_source_search"] == 2
        assert prog_a["queried_project_id"] == 2
        assert prog_a["grant_id"] == 3
        assert prog_a["earliest_fiscal_year"] == 2015

        # PROG_B: 1 unique search, 2 unique projects, 2 unique grants, min FY 2019
        prog_b = result[result["program.program_id"] == "PROG_B"].iloc[0]
        assert prog_b["api_source_search"] == 1
        assert prog_b["queried_project_id"] == 2
        assert prog_b["grant_id"] == 2
        assert prog_b["earliest_fiscal_year"] == 2019

    def test_fiscal_year_column_is_renamed(self, sample_grants, tmp_path):
        """Verify fiscal year column is renamed with 'earliest_' prefix."""
        output_file = tmp_path / "grantsStatsByProgram.csv"
        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_AGG_FUNCS_BY_COL = {
                "api_source_search": "nunique",
                "queried_project_id": "nunique",
                "grant_id": "nunique",
                "fiscal_year": "min",
            }
            mock_cfg.STAT_FISCALYEAR_COL = "fiscal_year"
            mock_cfg.STAT_GRANTS_BY_PROGRAM_FILENAME = str(output_file)

            get_grant_stats_by_program(sample_grants)

        result = pd.read_csv(output_file)
        assert "earliest_fiscal_year" in result.columns
        assert "fiscal_year" not in result.columns


# ============================================================
# get_shared_projects_by_program_pair
# ============================================================

class TestGetSharedProjectsByProgramPair:
    """Tests for identifying shared projects between program pairs."""

    def test_shared_project_detected(self, sample_grants, tmp_path):
        """R01CA222222 is shared between PROG_A and PROG_B."""
        output_file = tmp_path / "sharedProjectsByProgramPair.csv"
        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_CORE_PROJECT_COL = "queried_project_id"
            mock_cfg.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME = str(
                output_file)

            get_shared_projects_by_program_pair(sample_grants)

        result = pd.read_csv(output_file)
        assert not result.empty
        assert "program_1" in result.columns
        assert "program_2" in result.columns

        # Find the pair containing both PROG_A and PROG_B
        pair_rows = result[
            (result["program_1"].isin(["PROG_A", "PROG_B"]))
            & (result["program_2"].isin(["PROG_A", "PROG_B"]))
        ]
        assert len(pair_rows) == 1
        assert pair_rows.iloc[0]["unique_core_project_count"] == 1

    def test_no_shared_projects_returns_notice(self, single_program_grants,
                                               tmp_path):
        """When only one program exists, no pairs are possible."""
        output_file = tmp_path / "sharedProjectsByProgramPair.csv"
        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_CORE_PROJECT_COL = "queried_project_id"
            mock_cfg.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME = str(
                output_file)

            get_shared_projects_by_program_pair(single_program_grants)

        result = pd.read_csv(output_file)
        assert "NOTICE" in result.columns

    def test_output_csv_is_written(self, sample_grants, tmp_path):
        """Verify that the shared-projects report CSV is created."""
        output_file = tmp_path / "sharedProjectsByProgramPair.csv"
        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_CORE_PROJECT_COL = "queried_project_id"
            mock_cfg.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME = str(
                output_file)

            get_shared_projects_by_program_pair(sample_grants)

        assert output_file.exists()

    def test_multiple_shared_projects(self, tmp_path):
        """Two projects shared between the same pair increases count."""
        grants = pd.DataFrame({
            "program.program_id": ["A", "A", "B", "B"],
            "queried_project_id": ["P1", "P2", "P1", "P2"],
            "grant_id": ["G1", "G2", "G3", "G4"],
        })
        output_file = tmp_path / "shared.csv"
        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_CORE_PROJECT_COL = "queried_project_id"
            mock_cfg.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME = str(
                output_file)

            get_shared_projects_by_program_pair(grants)

        result = pd.read_csv(output_file)
        pair = result[
            (result["program_1"].isin(["A", "B"]))
            & (result["program_2"].isin(["A", "B"]))
        ]
        assert pair.iloc[0]["unique_core_project_count"] == 2


# ============================================================
# build_program_project_stats (orchestrator)
# ============================================================

class TestBuildProgramProjectStats:
    """Tests for the main orchestrator function."""

    def test_creates_reports_directory(self, sample_grants, tmp_path):
        """Verify the reports directory is created if missing."""
        reports_dir = tmp_path / "reports" / "new_version"
        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.REPORTS_GATHERED_DIR = str(reports_dir)
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_AGG_FUNCS_BY_COL = {
                "api_source_search": "nunique",
                "queried_project_id": "nunique",
                "grant_id": "nunique",
                "fiscal_year": "min",
            }
            mock_cfg.STAT_FISCALYEAR_COL = "fiscal_year"
            mock_cfg.STAT_CORE_PROJECT_COL = "queried_project_id"
            mock_cfg.STAT_GRANTS_BY_PROGRAM_FILENAME = str(
                reports_dir / "grantsStatsByProgram.csv")
            mock_cfg.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME = str(
                reports_dir / "sharedProjectsByProgramPair.csv")

            build_program_project_stats(sample_grants)

        assert reports_dir.exists()

    def test_both_reports_are_generated(self, sample_grants, tmp_path):
        """Verify both report files are produced by the orchestrator."""
        reports_dir = tmp_path / "reports"
        reports_dir.mkdir()
        grants_file = reports_dir / "grantsStatsByProgram.csv"
        shared_file = reports_dir / "sharedProjectsByProgramPair.csv"

        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.REPORTS_GATHERED_DIR = str(reports_dir)
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_AGG_FUNCS_BY_COL = {
                "api_source_search": "nunique",
                "queried_project_id": "nunique",
                "grant_id": "nunique",
                "fiscal_year": "min",
            }
            mock_cfg.STAT_FISCALYEAR_COL = "fiscal_year"
            mock_cfg.STAT_CORE_PROJECT_COL = "queried_project_id"
            mock_cfg.STAT_GRANTS_BY_PROGRAM_FILENAME = str(grants_file)
            mock_cfg.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME = str(
                shared_file)

            build_program_project_stats(sample_grants)

        assert grants_file.exists()
        assert shared_file.exists()

    def test_existing_reports_directory_is_reused(self, sample_grants,
                                                  tmp_path):
        """Verify no error when reports directory already exists."""
        reports_dir = tmp_path / "reports"
        reports_dir.mkdir()

        with patch("modules.build_program_project_stats.config") as mock_cfg:
            mock_cfg.REPORTS_GATHERED_DIR = str(reports_dir)
            mock_cfg.PROGRAM_ID_FIELDNAME = "program.program_id"
            mock_cfg.STAT_AGG_FUNCS_BY_COL = {
                "api_source_search": "nunique",
                "queried_project_id": "nunique",
                "grant_id": "nunique",
                "fiscal_year": "min",
            }
            mock_cfg.STAT_FISCALYEAR_COL = "fiscal_year"
            mock_cfg.STAT_CORE_PROJECT_COL = "queried_project_id"
            mock_cfg.STAT_GRANTS_BY_PROGRAM_FILENAME = str(
                reports_dir / "grantsStatsByProgram.csv")
            mock_cfg.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME = str(
                reports_dir / "sharedProjectsByProgramPair.csv")

            # Should not raise
            build_program_project_stats(sample_grants)
