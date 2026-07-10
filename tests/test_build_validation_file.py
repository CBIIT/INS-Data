"""
test_build_validation_file.py
2026-07-06 ZD

Pytest test suite for the build_validation_file.py module.
"""

import os

import pandas as pd
import pytest

from modules.build_validation_file import (
    get_single_node_counts,
    get_all_node_counts,
    get_downstream_node_records,
    get_detail_page_url,
    get_single_filter_value_results,
    build_single_program_output_counts,
    build_single_project_output_counts,
    build_multi_program_output_counts,
    save_dataframes_to_excel,
)


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def sample_programs():
    return pd.DataFrame({
        "type": ["program", "program"],
        "program_id": ["PROG1", "PROG2"],
        "program_name": ["Program One", "Program Two"],
        "program_acronym": ["P1", "P2"],
        "cancer_type": ["Breast Cancer", "Lung Cancer;Breast Cancer"],
        "focus_area": ["Treatment", "Prevention"],
    })


@pytest.fixture
def sample_projects():
    return pd.DataFrame({
        "type": ["project"] * 4,
        "project_id": ["PRJ_A", "PRJ_B", "PRJ_A", "PRJ_C"],
        "project_title": ["Title A", "Title B", "Title A", "Title C"],
        "program.program_id": ["PROG1", "PROG1", "PROG2", "PROG2"],
    })


@pytest.fixture
def sample_grants():
    return pd.DataFrame({
        "type": ["grant"] * 5,
        "grant_id": ["G1", "G2", "G3", "G4", "G5"],
        "project.project_id": ["PRJ_A", "PRJ_A", "PRJ_B", "PRJ_C", "PRJ_C"],
    })


@pytest.fixture
def sample_publications():
    return pd.DataFrame({
        "type": ["publication"] * 4,
        "pmid": [111, 222, 333, 444],
        "project.project_id": ["PRJ_A", "PRJ_B", "PRJ_C", "PRJ_C"],
    })


# ============================================================
# get_single_node_counts
# ============================================================

class TestGetSingleNodeCounts:
    """Tests for counting unique/total node IDs."""

    def test_returns_correct_counts(self, sample_projects):
        node_type, unique, total, multi_row = get_single_node_counts(
            sample_projects)

        assert node_type == "project"
        assert unique == 3   # PRJ_A, PRJ_B, PRJ_C
        assert total == 4    # 4 rows total
        assert multi_row == 1  # PRJ_A appears in 2 rows

    def test_all_unique_returns_zero_multi(self, sample_grants):
        _, _, _, multi_row = get_single_node_counts(sample_grants)
        assert multi_row == 0


# ============================================================
# get_all_node_counts
# ============================================================

class TestGetAllNodeCounts:
    """Tests for aggregating counts across multiple node types."""

    def test_produces_row_per_dataframe(
            self, sample_programs, sample_projects, sample_grants):
        result = get_all_node_counts(
            [sample_programs, sample_projects, sample_grants])

        assert len(result) == 3
        assert list(result["node_type"]) == [
            "program", "project", "grant"]


# ============================================================
# get_downstream_node_records
# ============================================================

class TestGetDownstreamNodeRecords:
    """Tests for filtering downstream records by link ID."""

    def test_single_value_filter(self, sample_projects):
        result = get_downstream_node_records(
            sample_projects, "program.program_id", "PROG1")
        assert len(result) == 2
        assert set(result["project_id"]) == {"PRJ_A", "PRJ_B"}

    def test_list_value_filter(self, sample_grants):
        result = get_downstream_node_records(
            sample_grants, "project.project_id", ["PRJ_A", "PRJ_C"])
        assert len(result) == 4  # G1, G2 for PRJ_A + G4, G5 for PRJ_C

    def test_no_match_returns_empty(self, sample_projects):
        result = get_downstream_node_records(
            sample_projects, "program.program_id", "NONEXISTENT")
        assert len(result) == 0


# ============================================================
# get_detail_page_url
# ============================================================

class TestGetDetailPageUrl:
    """Tests for building INS detail page URLs."""

    def test_dev_tier(self):
        url = get_detail_page_url("PROG1", "program", "-dev")
        assert url == "https://studycatalog-dev.cancer.gov/#/program/PROG1"

    def test_prod_tier(self):
        url = get_detail_page_url("PROG1", "program", "")
        assert url == "https://studycatalog.cancer.gov/#/program/PROG1"

    def test_project_url(self):
        url = get_detail_page_url("PRJ_A", "project", "-qa")
        assert url == "https://studycatalog-qa.cancer.gov/#/project/PRJ_A"


# ============================================================
# get_single_filter_value_results
# ============================================================

class TestGetSingleFilterValueResults:
    """Tests for filtering programs by facet values."""

    def test_exact_match_in_semicolon_list(self, sample_programs):
        result = get_single_filter_value_results(
            "cancer_type", "Breast Cancer",
            sample_programs, "program_id")
        # Both programs have Breast Cancer
        assert set(result) == {"PROG1", "PROG2"}

    def test_no_partial_match(self, sample_programs):
        """'Breast' alone should not match 'Breast Cancer'."""
        result = get_single_filter_value_results(
            "cancer_type", "Breast",
            sample_programs, "program_id")
        assert len(result) == 0

    def test_single_value_field(self, sample_programs):
        result = get_single_filter_value_results(
            "focus_area", "Treatment",
            sample_programs, "program_id")
        assert result == ["PROG1"]


# ============================================================
# build_single_program_output_counts
# ============================================================

class TestBuildSingleProgramOutputCounts:
    """Tests for program-level output count summaries."""

    def test_produces_row_per_program(
            self, sample_programs, sample_projects,
            sample_grants, sample_publications):
        result = build_single_program_output_counts(
            sample_programs, sample_projects,
            sample_grants, sample_publications)

        assert len(result) == 2
        assert set(result["program_id"]) == {"PROG1", "PROG2"}

    def test_counts_are_correct_for_prog1(
            self, sample_programs, sample_projects,
            sample_grants, sample_publications):
        result = build_single_program_output_counts(
            sample_programs, sample_projects,
            sample_grants, sample_publications)

        prog1 = result[result["program_id"] == "PROG1"].iloc[0]
        # PROG1 has projects PRJ_A, PRJ_B
        assert prog1["projects"] == 2
        # PRJ_A has G1,G2; PRJ_B has G3
        assert prog1["grants"] == 3
        # PRJ_A has 111; PRJ_B has 222
        assert prog1["publications"] == 2


# ============================================================
# build_single_project_output_counts
# ============================================================

class TestBuildSingleProjectOutputCounts:
    """Tests for project-level output count summaries."""

    def test_produces_row_per_unique_project(
            self, sample_projects, sample_grants, sample_publications):
        result = build_single_project_output_counts(
            sample_projects, sample_grants, sample_publications)

        # 3 unique project IDs
        assert len(result) == 3

    def test_counts_are_correct_for_prj_c(
            self, sample_projects, sample_grants, sample_publications):
        result = build_single_project_output_counts(
            sample_projects, sample_grants, sample_publications)

        prj_c = result[result["project_id"] == "PRJ_C"].iloc[0]
        assert prj_c["grants"] == 2    # G4, G5
        assert prj_c["publications"] == 2  # 333, 444


# ============================================================
# build_multi_program_output_counts
# ============================================================

class TestBuildMultiProgramOutputCounts:
    """Tests for combined counts across multiple programs."""

    def test_deduplicates_shared_projects(
            self, sample_projects, sample_grants, sample_publications):
        """PRJ_A is in both PROG1 and PROG2 — should be counted once."""
        result = build_multi_program_output_counts(
            ["PROG1", "PROG2"],
            sample_projects, sample_grants, sample_publications)

        # PRJ_A, PRJ_B, PRJ_C = 3 unique projects
        assert result["projects"] == 3

    def test_single_program_matches_individual(
            self, sample_projects, sample_grants, sample_publications):
        result = build_multi_program_output_counts(
            ["PROG1"],
            sample_projects, sample_grants, sample_publications)
        assert result["projects"] == 2  # PRJ_A, PRJ_B


# ============================================================
# save_dataframes_to_excel
# ============================================================

class TestSaveDataframesToExcel:
    """Tests for Excel output generation."""

    def test_creates_excel_with_correct_tabs(self, tmp_path):
        output_file = str(tmp_path / "test_output.xlsx")
        df_dict = {
            "sheet1": pd.DataFrame({"a": [1, 2], "b": [3, 4]}),
            "sheet2": pd.DataFrame({"x": [5, 6]}),
        }
        save_dataframes_to_excel(df_dict, output_file)

        assert os.path.exists(output_file)

        # Verify both sheets are present and readable
        loaded = pd.read_excel(output_file, sheet_name=None)
        assert set(loaded.keys()) == {"sheet1", "sheet2"}
        assert len(loaded["sheet1"]) == 2
        assert len(loaded["sheet2"]) == 2
