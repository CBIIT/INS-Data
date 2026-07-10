"""
test_gather_ctd2_data.py
2025-12-01 ZD

Pytest test suite for the `gather_ctd2_data.py` module.
"""

import uuid

import pandas as pd
import pytest
from unittest.mock import patch

from modules.gather_ctd2_data import (
    get_composite_uuid5,
    gather_ctd2_datasets,
    gather_ctd2_filedata,
    gather_ctd2_data,
)


# ============================================================
# Shared Fixtures
# ============================================================

@pytest.fixture
def datasets_df():
    """Sample CTD2 datasets DataFrame matching the input CSV schema."""
    return pd.DataFrame({
        'type': ['ctd2_dataset', 'ctd2_dataset'],
        'dataset_uuid': ['', ''],
        'dataset_source_repo': ['CTD² Network', 'CTD² Network'],
        'dataset_source_id': ['ctd2_001', 'ctd2_002'],
        'dataset_title': [
            'Computational Human High-grade Glioblastoma Multiform (GBM) '
            'Interactome - miRNA (Post-transcriptional) Layer',
            'Direct Reversal of Glucocorticoid Resistance by AKT Inhibition '
            'in Acute Lymphoblastic Leukemia (T-ALL)',
        ],
        'description': [
            'The Human High-Grade Glioma Interactome (HGi) contains a '
            'genome-wide complement of molecular interactions that are '
            'Glioblastoma Multiforme (GBM)-specific.',
            'The goal of this project is to identify key druggable '
            'regulators of glucocorticoid resistance in T-ALL.',
        ],
        'experimental_approaches': ['', ''],
        'download_file_links': [
            'High_Grade_GBM_Interactome.zip',
            'Master_Regulator_Analysis_T-ALL.zip',
        ],
        'institute': ['Columbia University', 'Columbia University'],
        'PI_name': ['Andrea Califano, Ph.D.', 'Andrea Califano, Ph.D.'],
        'POC_name': ['Prem Subramaniam', 'Prem Subramaniam'],
        'POC_email': [
            'ps2536@cumc.columbia.edu',
            'ps2536@cumc.columbia.edu',
        ],
        'dataset_pmid': ['22000015', '24291004'],
        'assay_method': ['', ''],
        'study_type': [
            'microRNA target predictions',
            'Microarray gene expression',
        ],
        'primary_disease': ['Brain Cancer', 'Leukemia'],
        'participant_count': ['', ''],
        'study_links': [
            '',
            'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32215',
        ],
        'related_genes': ['', ''],
    })


@pytest.fixture
def filedata_df():
    """Sample CTD2 file metadata DataFrame."""
    return pd.DataFrame({
        'type': ['file', 'file', 'file'],
        'file_id': ['', '', ''],
        'file_name': [
            'High_Grade_GBM_Interactome.zip',
            'Master_Regulator_Analysis_T-ALL.zip',
            'Nonexistent_File.zip',
        ],
        'file_type': ['ZIP', 'ZIP', 'ZIP'],
        'file_url': [
            'Columbia/High_Grade_GBM_Interactome.zip',
            'Columbia/Master_Regulator_Analysis_T-ALL.zip',
            'Columbia/Nonexistent_File.zip',
        ],
        'access_level': ['Open', 'Open', 'Open'],
    })


# ============================================================
# UUID5 Generation
# ============================================================

class TestGetCompositeUuid5:

    def test_duplicate_detection(self, datasets_df):
        df = pd.concat([datasets_df, datasets_df.iloc[[0]]],
                       ignore_index=True)
        with pytest.raises(ValueError, match="Duplicate UUID5 values found"):
            get_composite_uuid5(
                df,
                ['dataset_source_repo', 'dataset_title', 'description'],
                'dataset_uuid')

    def test_deterministic(self, filedata_df):
        fields = ['file_name', 'file_type', 'access_level']
        df1 = get_composite_uuid5(filedata_df.copy(), fields, 'file_id')
        df2 = get_composite_uuid5(filedata_df.copy(), fields, 'file_id')
        assert (df1['file_id'] == df2['file_id']).all()


# ============================================================
# gather_ctd2_datasets
# ============================================================

class TestGatherCtd2Datasets:

    @patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
    @patch('modules.gather_ctd2_data.pd.read_csv')
    @patch('modules.gather_ctd2_data.config')
    def test_produces_unique_uuids(self, mock_config, mock_read_csv,
                                    mock_to_csv, datasets_df, tmp_path):
        mock_config.CTD2_DATASET_INPUT_CSV = 'input.csv'
        mock_config.CTD2_DATASET_INTERMED_CSV = str(tmp_path / 'output.csv')
        mock_read_csv.return_value = datasets_df.copy()

        result = gather_ctd2_datasets()

        mock_to_csv.assert_called_once()
        assert 'dataset_uuid' in result.columns
        assert result['dataset_uuid'].nunique() == len(result)
        for uid in result['dataset_uuid']:
            uuid.UUID(str(uid))


# ============================================================
# gather_ctd2_filedata
# ============================================================

class TestGatherCtd2Filedata:

    @patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
    @patch('modules.gather_ctd2_data.pd.read_csv')
    @patch('modules.gather_ctd2_data.config')
    def test_produces_unique_file_ids(self, mock_config, mock_read_csv,
                                      mock_to_csv, filedata_df):
        mock_config.CTD2_FILE_INPUT_CSV = 'input.csv'
        mock_read_csv.return_value = filedata_df.copy()

        result = gather_ctd2_filedata()

        mock_to_csv.assert_not_called()
        assert 'file_id' in result.columns
        assert result['file_id'].nunique() == len(result)
        for fid in result['file_id']:
            uuid.UUID(str(fid))


# ============================================================
# Integration: gather_ctd2_data
# ============================================================

class TestGatherCtd2Data:

    @patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
    @patch('modules.gather_ctd2_data.os.makedirs')
    @patch('modules.gather_ctd2_data.config')
    @patch('modules.gather_ctd2_data.gather_ctd2_filedata')
    @patch('modules.gather_ctd2_data.gather_ctd2_datasets')
    def test_maps_files_to_datasets(self, mock_gather_datasets,
                                     mock_gather_filedata, mock_config,
                                     mock_makedirs, mock_to_csv,
                                     datasets_df, filedata_df):
        datasets = get_composite_uuid5(
            datasets_df.copy(),
            ['dataset_source_repo', 'dataset_title', 'description'],
            'dataset_uuid')
        files = get_composite_uuid5(
            filedata_df.copy(),
            ['file_name', 'file_type', 'access_level'],
            'file_id')
        mock_gather_datasets.return_value = datasets
        mock_gather_filedata.return_value = files
        mock_config.CTD2_FILE_INTERMED_CSV = 'output.csv'

        result = gather_ctd2_data()

        mock_to_csv.assert_called_once()
        assert isinstance(result, pd.DataFrame)
        assert 'dataset.dataset_uuid' in result.columns
        mapped = result[result['dataset.dataset_uuid'].notna()]
        unmapped = result[result['dataset.dataset_uuid'].isna()]
        assert set(mapped['file_name']) == {
            'High_Grade_GBM_Interactome.zip',
            'Master_Regulator_Analysis_T-ALL.zip',
        }
        assert set(unmapped['file_name']) == {'Nonexistent_File.zip'}

    @patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
    def test_empty_input(self, mock_to_csv, monkeypatch):
        cols_datasets = [
            'type', 'dataset_uuid', 'dataset_source_repo',
            'dataset_source_id', 'dataset_title', 'description',
            'experimental_approaches', 'download_file_links', 'institute',
            'PI_name', 'POC_name', 'POC_email', 'dataset_pmid',
            'assay_method', 'study_type', 'primary_disease',
            'participant_count', 'study_links', 'related_genes',
        ]
        cols_files = [
            'type', 'file_id', 'file_name', 'file_type', 'file_url',
            'access_level',
        ]
        empty_ds = pd.DataFrame(
            {c: pd.Series([], dtype='str') for c in cols_datasets})
        empty_fl = pd.DataFrame(
            {c: pd.Series([], dtype='str') for c in cols_files})
        monkeypatch.setattr(
            'modules.gather_ctd2_data.gather_ctd2_datasets', lambda: empty_ds)
        monkeypatch.setattr(
            'modules.gather_ctd2_data.gather_ctd2_filedata', lambda: empty_fl)

        result = gather_ctd2_data()

        mock_to_csv.assert_called_once()
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    @patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
    def test_missing_file_name_raises(self, mock_to_csv, monkeypatch):
        monkeypatch.setattr(
            'modules.gather_ctd2_data.gather_ctd2_datasets',
            lambda: pd.DataFrame({
                'type': ['ctd2_dataset'],
                'dataset_uuid': ['uuid1'],
                'dataset_source_repo': ['CTD² Network'],
                'dataset_source_id': ['ctd2_001'],
                'dataset_title': ['TitleA'],
                'description': ['DescA'],
                'experimental_approaches': [''],
                'download_file_links': ['High_Grade_GBM_Interactome.zip'],
                'institute': ['Columbia University'],
                'PI_name': ['Andrea Califano, Ph.D.'],
                'POC_name': ['Prem Subramaniam'],
                'POC_email': ['ps2536@cumc.columbia.edu'],
                'dataset_pmid': ['22000015'],
                'assay_method': [''],
                'study_type': ['microRNA target predictions'],
                'primary_disease': ['Brain Cancer'],
                'participant_count': [''],
                'study_links': [''],
                'related_genes': [''],
            }))
        monkeypatch.setattr(
            'modules.gather_ctd2_data.gather_ctd2_filedata',
            lambda: pd.DataFrame({
                'type': ['file'],
                'file_type': ['ZIP'],
                'access_level': ['Open'],
                'file_id': ['id1'],
                'file_url': ['Columbia/High_Grade_GBM_Interactome.zip'],
            }))

        with pytest.raises(KeyError):
            gather_ctd2_data()