"""
test_merge_curated_dbgap.py
2026-07-06 ZD

Pytest test suite for the merge_curated_dbgap.py module.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.merge_curated_dbgap import (
    load_dbgap_tsv,
    detect_title_changes,
    apply_title_decisions,
    merge_dbgap_datasets,
)


# ============================================================
# Fixtures
# ============================================================

COLUMNS = ["dataset_source_id", "dataset_title", "description",
           "dataset_source_repo", "dataset_uuid"]


@pytest.fixture
def old_curated_tsv(tmp_path):
    """Write an old curated TSV and return its path."""
    df = pd.DataFrame({
        "dataset_source_id": ["phs000001", "phs000002", "phs000003"],
        "dataset_title": ["Old Title 1", "Old Title 2", "Old Title 3"],
        "description": ["Desc 1", "Desc 2", "Desc 3"],
        "dataset_source_repo": ["dbGaP", "dbGaP", "dbGaP"],
        "dataset_uuid": ["uuid-1", "uuid-2", "uuid-3"],
    })
    path = tmp_path / "old_curated.tsv"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def new_gathered_tsv(tmp_path):
    """Write a new gathered TSV and return its path.
    - phs000001: shared, title changed
    - phs000002: shared, title unchanged
    - phs000004: new study
    - phs000003: removed (not in new)
    """
    df = pd.DataFrame({
        "dataset_source_id": ["phs000001", "phs000002", "phs000004"],
        "dataset_title": ["New Title 1", "Old Title 2", "New Study 4"],
        "description": ["Desc 1 new", "Desc 2", "Desc 4"],
        "dataset_source_repo": ["dbGaP", "dbGaP", "dbGaP"],
        "dataset_uuid": ["uuid-new-1", "uuid-new-2", "uuid-new-4"],
    })
    path = tmp_path / "new_gathered.tsv"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


# ============================================================
# load_dbgap_tsv
# ============================================================

class TestLoadDbgapTsv:
    """Tests for loading and validating dbGaP TSV files."""

    def test_loads_valid_tsv(self, old_curated_tsv):
        df = load_dbgap_tsv(old_curated_tsv)
        assert len(df) == 3
        assert "dataset_source_id" in df.columns

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_dbgap_tsv("/nonexistent/path.tsv")

    def test_missing_merge_key_raises(self, tmp_path):
        path = tmp_path / "bad.tsv"
        pd.DataFrame({"wrong_col": ["a"]}).to_csv(
            path, sep="\t", index=False)
        with pytest.raises(KeyError, match="dataset_source_id"):
            load_dbgap_tsv(str(path))


# ============================================================
# detect_title_changes
# ============================================================

class TestDetectTitleChanges:
    """Tests for detecting title differences between old and new datasets."""

    def test_detects_changed_titles(self):
        old = pd.DataFrame({
            "dataset_source_id": ["phs000001", "phs000002"],
            "dataset_title": ["Old Title", "Same Title"],
        })
        new = pd.DataFrame({
            "dataset_source_id": ["phs000001", "phs000002"],
            "dataset_title": ["New Title", "Same Title"],
        })
        result = detect_title_changes(old, new)
        assert len(result) == 1
        assert result.iloc[0]["old_title"] == "Old Title"
        assert result.iloc[0]["new_title"] == "New Title"

    def test_no_changes_returns_empty(self):
        old = pd.DataFrame({
            "dataset_source_id": ["phs000001"],
            "dataset_title": ["Same"],
        })
        new = old.copy()
        result = detect_title_changes(old, new)
        assert len(result) == 0


# ============================================================
# apply_title_decisions
# ============================================================

class TestApplyTitleDecisions:
    """Tests for applying user title change decisions."""

    def test_applies_accepted_changes(self, tmp_path):
        merged = pd.DataFrame({
            "dataset_source_id": ["phs000001", "phs000002"],
            "dataset_title": ["Old Title 1", "Old Title 2"],
        })
        review = pd.DataFrame({
            "dataset_source_id": ["phs000001", "phs000002"],
            "old_title": ["Old Title 1", "Old Title 2"],
            "new_title": ["New Title 1", "New Title 2"],
            "use_new_title": ["yes", "no"],
        })
        review_path = tmp_path / "review.csv"
        review.to_csv(review_path, index=False)

        result = apply_title_decisions(merged, str(review_path))

        assert result.iloc[0]["dataset_title"] == "New Title 1"
        assert result.iloc[1]["dataset_title"] == "Old Title 2"

    def test_no_acceptances_returns_unchanged(self, tmp_path):
        merged = pd.DataFrame({
            "dataset_source_id": ["phs000001"],
            "dataset_title": ["Old Title"],
        })
        review = pd.DataFrame({
            "dataset_source_id": ["phs000001"],
            "old_title": ["Old Title"],
            "new_title": ["New Title"],
            "use_new_title": ["no"],
        })
        review_path = tmp_path / "review.csv"
        review.to_csv(review_path, index=False)

        result = apply_title_decisions(merged, str(review_path))
        assert result.iloc[0]["dataset_title"] == "Old Title"


# ============================================================
# merge_dbgap_datasets (integration)
# ============================================================

class TestMergeDbgapDatasets:
    """Integration tests for the full merge workflow."""

    def test_merge_preserves_old_and_appends_new(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        output = str(tmp_path / "merged.tsv")
        review = str(tmp_path / "review.csv")

        result = merge_dbgap_datasets(
            old_curated_tsv, new_gathered_tsv, output, review)

        # Old: phs000001, phs000002, phs000003
        # New only: phs000004
        # Total: 4
        assert len(result) == 4

        ids = set(result["dataset_source_id"])
        assert ids == {"phs000001", "phs000002", "phs000003", "phs000004"}

    def test_old_rows_take_precedence(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        """For shared studies, old curated values should be kept."""
        output = str(tmp_path / "merged.tsv")
        review = str(tmp_path / "review.csv")

        result = merge_dbgap_datasets(
            old_curated_tsv, new_gathered_tsv, output, review)

        phs1 = result[
            result["dataset_source_id"] == "phs000001"].iloc[0]
        # Old title should be kept (no review CSV applied yet)
        assert phs1["dataset_title"] == "Old Title 1"

    def test_new_only_rows_marked_as_new(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        output = str(tmp_path / "merged.tsv")
        review = str(tmp_path / "review.csv")

        result = merge_dbgap_datasets(
            old_curated_tsv, new_gathered_tsv, output, review)

        phs4 = result[
            result["dataset_source_id"] == "phs000004"].iloc[0]
        assert phs4["curation_status"] == "new_study"

    def test_old_only_rows_retained(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        """phs000003 is only in old — should be retained."""
        output = str(tmp_path / "merged.tsv")
        review = str(tmp_path / "review.csv")

        result = merge_dbgap_datasets(
            old_curated_tsv, new_gathered_tsv, output, review)

        assert "phs000003" in result["dataset_source_id"].values

    def test_exports_title_review_csv(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        """Title changes should produce a review CSV on first run."""
        output = str(tmp_path / "merged.tsv")
        review = str(tmp_path / "review.csv")

        merge_dbgap_datasets(
            old_curated_tsv, new_gathered_tsv, output, review)

        assert os.path.exists(review)
        review_df = pd.read_csv(review)
        # phs000001 has a title change
        assert len(review_df) == 1
        assert review_df.iloc[0]["dataset_source_id"] == "phs000001"

    def test_output_file_created(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        output = str(tmp_path / "merged.tsv")
        review = str(tmp_path / "review.csv")

        merge_dbgap_datasets(
            old_curated_tsv, new_gathered_tsv, output, review)

        assert os.path.exists(output)
        loaded = pd.read_csv(output, sep="\t")
        assert len(loaded) == 4
