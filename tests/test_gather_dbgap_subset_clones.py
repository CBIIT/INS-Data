"""
test_gather_dbgap_subset_clones.py
2026-07-06 ZD

Pytest test suite for the gather_dbgap_subset_clones.py module.
"""

import os
import sys
import uuid

import pandas as pd
import pytest
from unittest.mock import patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.gather_dbgap_subset_clones import (
    load_subset_phs,
    generate_clone_uuid,
    build_subset_clones,
)


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def sample_dbgap_df():
    """Minimal curated clean dbGaP DataFrame."""
    return pd.DataFrame({
        "dataset_source_id": ["phs000001", "phs000002", "phs000003"],
        "dataset_title": ["Study A", "Study B", "Study C"],
        "dataset_source_repo": ["dbGaP", "dbGaP", "dbGaP"],
        "dataset_uuid": ["uuid-1", "uuid-2", "uuid-3"],
        "dataset_storage_distribution": ["CRDC", "GDC", ""],
    })


# ============================================================
# generate_clone_uuid
# ============================================================

class TestGenerateCloneUuid:
    """Tests for deterministic UUID generation for clones."""

    def test_produces_valid_uuid(self):
        result = generate_clone_uuid("CRDC", "phs000001")
        # Should be a valid UUID string
        uuid.UUID(result)

    def test_deterministic(self):
        """Same inputs should always produce the same UUID."""
        a = generate_clone_uuid("CRDC", "phs000001")
        b = generate_clone_uuid("CRDC", "phs000001")
        assert a == b

    def test_different_repos_produce_different_uuids(self):
        crdc = generate_clone_uuid("CRDC", "phs000001")
        gdc = generate_clone_uuid("GDC", "phs000001")
        assert crdc != gdc

    def test_different_phs_produce_different_uuids(self):
        a = generate_clone_uuid("CRDC", "phs000001")
        b = generate_clone_uuid("CRDC", "phs000002")
        assert a != b


# ============================================================
# load_subset_phs
# ============================================================

class TestLoadSubsetPhs:
    """Tests for loading subset CSV files."""

    def test_loads_valid_csv(self, tmp_path):
        csv_path = tmp_path / "CRDC_subset_2026-03-09.csv"
        pd.DataFrame({
            "accession": ["phs000001.v1.p1", "phs000002.v3.p2"]
        }).to_csv(csv_path, index=False)

        with patch("modules.gather_dbgap_subset_clones.config") as mock:
            mock.DBGAP_SUBSET_INPUT_DIR = str(tmp_path)
            mock.DBGAP_CSV_VERSION = "2026-03-09"
            result = load_subset_phs("CRDC")

        assert result == {"phs000001", "phs000002"}

    def test_strips_versioning(self, tmp_path):
        csv_path = tmp_path / "GDC_subset_2026-03-09.csv"
        pd.DataFrame({
            "accession": ["phs999999.v9.p3"]
        }).to_csv(csv_path, index=False)

        with patch("modules.gather_dbgap_subset_clones.config") as mock:
            mock.DBGAP_SUBSET_INPUT_DIR = str(tmp_path)
            mock.DBGAP_CSV_VERSION = "2026-03-09"
            result = load_subset_phs("GDC")

        assert result == {"phs999999"}

    def test_missing_file_returns_empty(self, tmp_path):
        with patch("modules.gather_dbgap_subset_clones.config") as mock:
            mock.DBGAP_SUBSET_INPUT_DIR = str(tmp_path)
            mock.DBGAP_CSV_VERSION = "2026-03-09"
            result = load_subset_phs("NONEXISTENT")

        assert result == set()


# ============================================================
# build_subset_clones
# ============================================================

class TestBuildSubsetClones:
    """Tests for the full clone-building logic."""

    def test_clones_matching_rows(self, sample_dbgap_df, tmp_path):
        """Rows matching the subset CSV should be cloned with new repo+UUID."""
        # Create a CRDC subset CSV with phs000001
        csv_path = tmp_path / "CRDC_subset_2026-03-09.csv"
        pd.DataFrame({
            "accession": ["phs000001.v1.p1"]
        }).to_csv(csv_path, index=False)

        with patch("modules.gather_dbgap_subset_clones.config") as mock:
            mock.DBGAP_SUBSET_INPUT_DIR = str(tmp_path)
            mock.DBGAP_CSV_VERSION = "2026-03-09"
            mock.DBGAP_STORAGE_DISTRIBUTION_MAP = {"CRDC": "CRDC"}
            mock.DBGAP_SUBSET_REPO_MAP = {"CRDC": "CRDC"}

            result = build_subset_clones(sample_dbgap_df)

        assert len(result) == 1
        assert result.iloc[0]["dataset_source_repo"] == "CRDC"
        assert result.iloc[0]["dataset_source_id"] == "phs000001"
        # UUID should be different from original
        assert result.iloc[0]["dataset_uuid"] != "uuid-1"

    def test_preserves_source_id_and_other_fields(
            self, sample_dbgap_df, tmp_path):
        """Clone should keep dataset_source_id and other fields unchanged."""
        csv_path = tmp_path / "GDC_subset_2026-03-09.csv"
        pd.DataFrame({
            "accession": ["phs000002"]
        }).to_csv(csv_path, index=False)

        with patch("modules.gather_dbgap_subset_clones.config") as mock:
            mock.DBGAP_SUBSET_INPUT_DIR = str(tmp_path)
            mock.DBGAP_CSV_VERSION = "2026-03-09"
            mock.DBGAP_STORAGE_DISTRIBUTION_MAP = {"GDC": "GDC"}
            mock.DBGAP_SUBSET_REPO_MAP = {"GDC": "GDC"}

            result = build_subset_clones(sample_dbgap_df)

        row = result.iloc[0]
        assert row["dataset_source_id"] == "phs000002"
        assert row["dataset_title"] == "Study B"

    def test_multiple_subsets_same_repo_deduplicates(
            self, sample_dbgap_df, tmp_path):
        """CRDC and CTDC both mapping to 'CRDC' should deduplicate."""
        # Both subsets reference phs000001
        for key in ("CRDC", "CTDC"):
            csv_path = tmp_path / f"{key}_subset_2026-03-09.csv"
            pd.DataFrame({
                "accession": ["phs000001"]
            }).to_csv(csv_path, index=False)

        with patch("modules.gather_dbgap_subset_clones.config") as mock:
            mock.DBGAP_SUBSET_INPUT_DIR = str(tmp_path)
            mock.DBGAP_CSV_VERSION = "2026-03-09"
            mock.DBGAP_STORAGE_DISTRIBUTION_MAP = {
                "CRDC": "CRDC", "CTDC": "CRDC"}
            mock.DBGAP_SUBSET_REPO_MAP = {
                "CRDC": "CRDC", "CTDC": "CRDC"}

            result = build_subset_clones(sample_dbgap_df)

        # Should only produce one clone for phs000001 under CRDC
        assert len(result) == 1

    def test_no_matching_rows_returns_empty(
            self, sample_dbgap_df, tmp_path):
        csv_path = tmp_path / "CRDC_subset_2026-03-09.csv"
        pd.DataFrame({
            "accession": ["phs999999"]
        }).to_csv(csv_path, index=False)

        with patch("modules.gather_dbgap_subset_clones.config") as mock:
            mock.DBGAP_SUBSET_INPUT_DIR = str(tmp_path)
            mock.DBGAP_CSV_VERSION = "2026-03-09"
            mock.DBGAP_STORAGE_DISTRIBUTION_MAP = {"CRDC": "CRDC"}
            mock.DBGAP_SUBSET_REPO_MAP = {"CRDC": "CRDC"}

            result = build_subset_clones(sample_dbgap_df)

        assert len(result) == 0
