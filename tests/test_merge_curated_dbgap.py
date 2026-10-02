"""Tests for the version-aware dbGaP curation merge."""

import os

import pandas as pd
import pytest

from modules.merge_curated_dbgap import load_dbgap_tsv, merge_dbgap_datasets


@pytest.fixture
def old_curated_tsv(tmp_path):
    df = pd.DataFrame({
        "dataset_source_id": ["phs000001", "phs000002", "phs000003"],
        "dataset_source_accession": [
            "phs000001.v1.p1", "phs000002.v1.p1", "phs000003.v1.p1"],
        "dataset_title": ["Curated Title 1", "Curated Title 2",
                          "Curated Title 3"],
        "description": ["Curated 1", "Curated 2", "Curated 3"],
        "dataset_source_repo": ["dbGaP", "dbGaP", "dbGaP"],
        "dataset_uuid": ["uuid-1", "uuid-2", "uuid-3"],
    })
    path = tmp_path / "old_curated.tsv"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def new_gathered_tsv(tmp_path):
    df = pd.DataFrame({
        "dataset_source_id": ["phs000001", "phs000002", "phs000004"],
        "dataset_source_accession": [
            "phs000001.v2.p1", "phs000002.v1.p1", "phs000004.v1.p1"],
        "dataset_title": ["New Title 1", "Source Title 2", "New Study 4"],
        "description": ["New 1", "Source 2", "New 4"],
        "dataset_source_repo": ["dbGaP", "dbGaP", "dbGaP"],
        "dataset_uuid": ["uuid-new-1", "uuid-new-2", "uuid-new-4"],
    })
    path = tmp_path / "new_gathered.tsv"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


class TestLoadDbgapTsv:
    def test_loads_valid_tsv(self, old_curated_tsv):
        assert len(load_dbgap_tsv(old_curated_tsv)) == 3

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_dbgap_tsv("/nonexistent/path.tsv")

    @pytest.mark.parametrize(
        "missing_column", ["dataset_source_id", "dataset_source_accession"])
    def test_missing_required_column_raises(
            self, old_curated_tsv, tmp_path, missing_column):
        df = pd.read_csv(old_curated_tsv, sep="\t").drop(
            columns=missing_column)
        path = tmp_path / "bad.tsv"
        df.to_csv(path, sep="\t", index=False)

        with pytest.raises(KeyError, match=missing_column):
            load_dbgap_tsv(str(path))

    def test_duplicate_short_accession_raises(self, old_curated_tsv, tmp_path):
        df = pd.read_csv(old_curated_tsv, sep="\t")
        df = pd.concat([df, df.iloc[[0]]], ignore_index=True)
        path = tmp_path / "duplicates.tsv"
        df.to_csv(path, sep="\t", index=False)

        with pytest.raises(ValueError, match="phs000001"):
            load_dbgap_tsv(str(path))


class TestMergeDbgapDatasets:
    def merge(self, old_path, new_path, tmp_path):
        output = str(tmp_path / "merged.tsv")
        result = merge_dbgap_datasets(old_path, new_path, output)
        return result, output

    def test_classifies_and_orders_rows(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        result, _ = self.merge(
            old_curated_tsv, new_gathered_tsv, tmp_path)

        assert list(result["dataset_source_id"]) == [
            "phs000001", "phs000004", "phs000003", "phs000002"]
        assert list(result["curation_status"]) == [
            "updated", "new", "old_only", "unchanged"]

    def test_updated_row_uses_new_gathered_data(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        result, _ = self.merge(
            old_curated_tsv, new_gathered_tsv, tmp_path)
        row = result[result["dataset_source_id"] == "phs000001"].iloc[0]

        assert row["dataset_title"] == "New Title 1"
        assert row["description"] == "New 1"
        assert row["previous_source_accession"] == "phs000001.v1.p1"
        assert row["current_source_accession"] == "phs000001.v2.p1"

    def test_unchanged_row_preserves_curated_data(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        result, _ = self.merge(
            old_curated_tsv, new_gathered_tsv, tmp_path)
        row = result[result["dataset_source_id"] == "phs000002"].iloc[0]

        assert row["dataset_title"] == "Curated Title 2"
        assert row["description"] == "Curated 2"
        assert row["previous_source_accession"] == "phs000002.v1.p1"
        assert row["current_source_accession"] == "phs000002.v1.p1"

    def test_new_row_uses_new_data(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        result, _ = self.merge(
            old_curated_tsv, new_gathered_tsv, tmp_path)
        row = result[result["dataset_source_id"] == "phs000004"].iloc[0]

        assert row["dataset_title"] == "New Study 4"
        assert row["previous_source_accession"] == ""
        assert row["current_source_accession"] == "phs000004.v1.p1"

    def test_old_only_row_preserves_curated_data(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        result, _ = self.merge(
            old_curated_tsv, new_gathered_tsv, tmp_path)
        row = result[result["dataset_source_id"] == "phs000003"].iloc[0]

        assert row["dataset_title"] == "Curated Title 3"
        assert row["previous_source_accession"] == "phs000003.v1.p1"
        assert row["current_source_accession"] == ""

    def test_output_file_created(
            self, old_curated_tsv, new_gathered_tsv, tmp_path):
        result, output = self.merge(
            old_curated_tsv, new_gathered_tsv, tmp_path)

        assert os.path.exists(output)
        loaded = pd.read_csv(output, sep="\t", keep_default_na=False)
        assert len(loaded) == len(result) == 4