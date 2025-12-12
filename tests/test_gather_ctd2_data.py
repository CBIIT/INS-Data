"""
test_gather_ctd2_data.py
2025-12-01 ZD

Pytest test suite for the `gather_ctd2_data.py` module.
"""

import os
import sys
import pandas as pd
import pytest
from unittest.mock import patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.gather_ctd2_data import (
    get_composite_uuid5,
    gather_ctd2_datasets,
    gather_ctd2_filedata,
    gather_ctd2_data,
)

@pytest.fixture
def datasets_df():
    return pd.DataFrame({
        'type': ['ctd2_dataset', 'ctd2_dataset'],
        'dataset_uuid': ['', ''],
        'dataset_source_repo': ['CTD² Network', 'CTD² Network'],
        'dataset_source_id': ['ctd2_001', 'ctd2_002'],
        'dataset_title': [
            'Computational Human High-grade Glioblastoma Multiform (GBM) Interactome - miRNA (Post-transcriptional) Layer',
            'Direct Reversal of Glucocorticoid Resistance by AKT Inhibition in Acute Lymphoblastic Leukemia (T-ALL)'
        ],
        'description': [
            'The Human High-Grade Glioma Interactome (HGi) contains a genome-wide complement of molecular interactions that are Glioblastoma Multiforme (GBM)-specific.',
            'The goal of this project is to identify key druggable regulators of glucocorticoid resistance in T-ALL.'
        ],
        'experimental_approaches': ['', ''],
        'download_file_links': ['High_Grade_GBM_Interactome.zip', 'Master_Regulator_Analysis_T-ALL.zip'],
        'institute': ['Columbia University', 'Columbia University'],
        'PI_name': ['Andrea Califano, Ph.D.', 'Andrea Califano, Ph.D.'],
        'POC_name': ['Prem Subramaniam', 'Prem Subramaniam'],
        'POC_email': ['ps2536@cumc.columbia.edu', 'ps2536@cumc.columbia.edu'],
        'dataset_pmid': ['22000015', '24291004'],
        'assay_method': ['', ''],
        'study_type': ['microRNA target predictions', 'Microarray gene expression'],
        'primary_disease': ['Brain Cancer', 'Leukemia'],
        'participant_count': ['', ''],
        'study_links': ['', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32215'],
        'related_genes': ['', '']
    })

@pytest.fixture
def filedata_df():
    return pd.DataFrame({
        'type': ['file', 'file', 'file'],
        'file_id': ['', '', ''],
        'file_name': [
            'High_Grade_GBM_Interactome.zip',
            'Master_Regulator_Analysis_T-ALL.zip',
            'Nonexistent_File.zip'
        ],
        'file_type': ['ZIP', 'ZIP', 'ZIP'],
        'file_url': [
            'Columbia/High_Grade_GBM_Interactome.zip',
            'Columbia/Master_Regulator_Analysis_T-ALL.zip',
            'Columbia/Nonexistent_File.zip'
        ],
        'access_level': ['Open', 'Open', 'Open']
    })

@patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
@patch('modules.gather_ctd2_data.pd.read_csv')
@patch('modules.gather_ctd2_data.config')
def test_gather_ctd2_datasets_basic(mock_config, mock_read_csv, mock_to_csv, datasets_df, tmp_path):
    mock_config.CTD2_DATASET_INPUT_CSV = 'input.csv'
    mock_config.CTD2_DATASET_INTERMED_CSV = str(tmp_path / 'output.csv')
    mock_read_csv.return_value = datasets_df.copy()
    result = gather_ctd2_datasets()
    mock_to_csv.assert_called_once()
    assert 'dataset_uuid' in result.columns
    assert result['dataset_uuid'].nunique() == len(result)

@patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
@patch('modules.gather_ctd2_data.pd.read_csv')
@patch('modules.gather_ctd2_data.config')
def test_gather_ctd2_filedata_basic(mock_config, mock_read_csv, mock_to_csv, filedata_df):
    mock_config.CTD2_FILE_INPUT_CSV = 'input.csv'
    mock_read_csv.return_value = filedata_df.copy()
    result = gather_ctd2_filedata()
    mock_to_csv.assert_not_called()  # gather_ctd2_filedata does not save directly
    assert 'file_id' in result.columns
    assert result['file_id'].nunique() == len(result)

@patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
@patch('modules.gather_ctd2_data.os.makedirs')
@patch('modules.gather_ctd2_data.config')
@patch('modules.gather_ctd2_data.gather_ctd2_filedata')
@patch('modules.gather_ctd2_data.gather_ctd2_datasets')
def test_gather_ctd2_data_mapping(
    mock_gather_datasets, mock_gather_filedata, mock_config, mock_makedirs, mock_to_csv,
    datasets_df, filedata_df
):
    datasets = get_composite_uuid5(
        datasets_df.copy(),
        ['dataset_source_repo', 'dataset_title', 'description'],
        'dataset_uuid'
    )
    files = get_composite_uuid5(
        filedata_df.copy(),
        ['file_name', 'file_type', 'access_level'],
        'file_id'
    )
    mock_gather_datasets.return_value = datasets
    mock_gather_filedata.return_value = files
    mock_config.CTD2_FILE_INTERMED_CSV = 'output.csv'
    result = gather_ctd2_data()
    mock_to_csv.assert_called_once()
    assert isinstance(result, pd.DataFrame)
    assert 'dataset.dataset_uuid' in result.columns
    mapped = result[result['dataset.dataset_uuid'].notna()]
    unmapped = result[result['dataset.dataset_uuid'].isna()]
    assert set(mapped['file_name']) == {'High_Grade_GBM_Interactome.zip', 'Master_Regulator_Analysis_T-ALL.zip'}
    assert set(unmapped['file_name']) == {'Nonexistent_File.zip'}

@patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
def test_get_composite_uuid5_duplicate_detection(mock_to_csv, datasets_df):
    df = pd.concat([datasets_df, datasets_df.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError):
        get_composite_uuid5(df, ['dataset_source_repo', 'dataset_title', 'description'], 'dataset_uuid')

@patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
def test_get_composite_uuid5_deterministic(mock_to_csv, filedata_df):
    df1 = get_composite_uuid5(filedata_df.copy(), ['file_name', 'file_type', 'access_level'], 'file_id')
    df2 = get_composite_uuid5(filedata_df.copy(), ['file_name', 'file_type', 'access_level'], 'file_id')
    assert (df1['file_id'] == df2['file_id']).all()

@patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
def test_gather_ctd2_data_empty(mock_to_csv, monkeypatch):
    empty_datasets = pd.DataFrame({
        'type': pd.Series([], dtype='str'),
        'dataset_uuid': pd.Series([], dtype='str'),
        'dataset_source_repo': pd.Series([], dtype='str'),
        'dataset_source_id': pd.Series([], dtype='str'),
        'dataset_title': pd.Series([], dtype='str'),
        'description': pd.Series([], dtype='str'),
        'experimental_approaches': pd.Series([], dtype='str'),
        'download_file_links': pd.Series([], dtype='str'),
        'institute': pd.Series([], dtype='str'),
        'PI_name': pd.Series([], dtype='str'),
        'POC_name': pd.Series([], dtype='str'),
        'POC_email': pd.Series([], dtype='str'),
        'dataset_pmid': pd.Series([], dtype='str'),
        'assay_method': pd.Series([], dtype='str'),
        'study_type': pd.Series([], dtype='str'),
        'primary_disease': pd.Series([], dtype='str'),
        'participant_count': pd.Series([], dtype='str'),
        'study_links': pd.Series([], dtype='str'),
        'related_genes': pd.Series([], dtype='str')
    })
    empty_filedata = pd.DataFrame({
        'type': pd.Series([], dtype='str'),
        'file_id': pd.Series([], dtype='str'),
        'file_name': pd.Series([], dtype='str'),
        'file_type': pd.Series([], dtype='str'),
        'file_url': pd.Series([], dtype='str'),
        'access_level': pd.Series([], dtype='str')
    })
    monkeypatch.setattr('modules.gather_ctd2_data.gather_ctd2_datasets', lambda: empty_datasets)
    monkeypatch.setattr('modules.gather_ctd2_data.gather_ctd2_filedata', lambda: empty_filedata)
    result = gather_ctd2_data()
    mock_to_csv.assert_called_once()
    assert isinstance(result, pd.DataFrame)
    assert result.empty

@patch('modules.gather_ctd2_data.pd.DataFrame.to_csv')
def test_gather_ctd2_data_missing_file_name(mock_to_csv, monkeypatch):
    monkeypatch.setattr('modules.gather_ctd2_data.gather_ctd2_datasets', lambda: pd.DataFrame({
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
        'related_genes': ['']
    }))
    monkeypatch.setattr('modules.gather_ctd2_data.gather_ctd2_filedata', lambda: pd.DataFrame({
        'type': ['file'],
        'file_type': ['ZIP'],
        'access_level': ['Open'],
        'file_id': ['id1'],
        'file_url': ['Columbia/High_Grade_GBM_Interactome.zip']
    }))
    with pytest.raises(KeyError):
        gather_ctd2_data()