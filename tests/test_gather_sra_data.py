"""
test_gather_sra_data.py
2026-02-13 ZD

Pytest test suite for the `gather_sra_data.py` module.
"""

import os
import sys
import json
import pandas as pd
import pytest
import uuid
from unittest.mock import patch, MagicMock
from io import BytesIO


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.gather_sra_data import (
    get_composite_uuid5,
    fetch_sra_ids,
    get_sra_ids_for_pubmed_ids,
    get_srp_ids_for_sra_id,
    get_srp_ids_for_sra_ids,
    merge_semicolon_fields,
    aggregate_batch_mappings,
    aggregate_batch_srp_data,
    get_max_batch_number,
    load_processed_pmids_from_batches,
    save_failed_pmids_report,
    load_failed_pmids_report,
    load_all_batch_files,
    _create_blank_metadata_dict,
    create_sra_mapping_dataframe,
    create_srp_centric_dataframe,
    fetch_sra_metadata_for_srp,
    enrich_srp_data_with_metadata,
    get_dataset_doc_from_project,
    create_srp_to_sra_mapping,
    gather_srp_metadata,
)


# ============================================================================
# Test Fixtures
# ============================================================================

@pytest.fixture
def sample_dataframe():
    """Create a sample DataFrame for testing UUID generation."""
    return pd.DataFrame({
        'dataset_source_id': ['SRP123', 'SRP456', 'SRP789'],
        'dataset_title': ['Study A', 'Study B', 'Study C'],
        'dataset_pmid': ['12345', '67890', '11111']
    })


@pytest.fixture
def sample_mapping_data():
    """Sample PMID to SRA ID mapping data."""
    return {
        '12345': ['SRX123', 'SRX456'],
        '67890': ['SRX789'],
        '11111': []
    }


@pytest.fixture
def sample_batch_dir(tmp_path):
    """Create a temporary batch directory with sample files.
    
    Uses uppercase column names (PMID, SRA_IDs, SRP_ERP_IDs) to match
    the column conventions used by create_sra_mapping_dataframe() and
    load_processed_pmids_from_batches().
    """
    batch_dir = tmp_path / "sra_batches"
    batch_dir.mkdir()
    
    # Create sample batch files with correct column names
    batch_1_mapping = pd.DataFrame({
        'PMID': ['12345', '67890'],
        'SRA_IDs': ['SRX123', 'SRX789'],
        'SRP_ERP_IDs': ['SRP001', 'SRP002']
    })
    batch_1_mapping.to_csv(batch_dir / "batch_0001_mapping.csv", index=False)
    
    batch_2_mapping = pd.DataFrame({
        'PMID': ['11111', '22222'],
        'SRA_IDs': ['SRX999', 'SRX888'],
        'SRP_ERP_IDs': ['SRP003', 'SRP004']
    })
    batch_2_mapping.to_csv(batch_dir / "batch_0002_mapping.csv", index=False)
    
    return str(batch_dir)


@pytest.fixture
def mock_entrez_response():
    """Mock successful Entrez elink response."""
    return [{
        'LinkSetDb': [{
            'Link': [
                {'Id': 'SRX123456'},
                {'Id': 'SRX789012'}
            ]
        }]
    }]


@pytest.fixture
def mock_xml_metadata():
    """Mock XML response for SRA metadata."""
    return """<?xml version="1.0" encoding="UTF-8"?>
    <EXPERIMENT_PACKAGE_SET>
        <EXPERIMENT_PACKAGE>
            <STUDY accession="SRP123456">
                <DESCRIPTOR>
                    <STUDY_TITLE>Test Study Title</STUDY_TITLE>
                    <STUDY_ABSTRACT>This is a test study abstract.</STUDY_ABSTRACT>
                    <STUDY_DESCRIPTION>Test study description with more details.</STUDY_DESCRIPTION>
                </DESCRIPTOR>
            </STUDY>
            <EXPERIMENT accession="SRX123456">
                <TITLE>Experiment Title</TITLE>
                <DESIGN>
                    <DESIGN_DESCRIPTION>Test design description</DESIGN_DESCRIPTION>
                    <SAMPLE_DESCRIPTOR accession="SRS123456">
                        <POOL>
                            <MEMBER accession="SAMN123456"/>
                        </POOL>
                    </SAMPLE_DESCRIPTOR>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
                        <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
                        <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
                        <LIBRARY_LAYOUT>
                            <PAIRED/>
                        </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                </DESIGN>
                <PLATFORM>
                    <ILLUMINA>
                        <INSTRUMENT_MODEL>Illumina HiSeq 2000</INSTRUMENT_MODEL>
                    </ILLUMINA>
                </PLATFORM>
            </EXPERIMENT>
            <SUBMISSION accession="SRA123456" center_name="TEST_CENTER">
                <CONTACTS>
                    <CONTACT name="John Doe"/>
                </CONTACTS>
            </SUBMISSION>
            <Organization type="center">
                <Name>Test Organization</Name>
            </Organization>
            <RUN_SET>
                <RUN accession="SRR123456" total_bases="1000000000" total_spots="10000000" published="2023-01-15">
                    <EXPERIMENT_REF accession="SRX123456"/>
                </RUN>
            </RUN_SET>
        </EXPERIMENT_PACKAGE>
    </EXPERIMENT_PACKAGE_SET>"""


# ============================================================================
# Test UUID5 Generation
# ============================================================================

def test_get_composite_uuid5_basic(sample_dataframe):
    """Test basic UUID5 generation from composite fields."""
    result = get_composite_uuid5(
        sample_dataframe, 
        fields=['dataset_source_id'],
        uuid_col='dataset_uuid'
    )
    
    # Check that UUID column was created
    assert 'dataset_uuid' in result.columns
    
    # Check that all UUIDs are valid UUID objects
    assert all(isinstance(uid, uuid.UUID) for uid in result['dataset_uuid'])
    
    # Check that UUIDs are deterministic (same input = same UUID)
    result2 = get_composite_uuid5(
        sample_dataframe.copy(), 
        fields=['dataset_source_id'],
        uuid_col='dataset_uuid'
    )
    assert result['dataset_uuid'].equals(result2['dataset_uuid'])


def test_get_composite_uuid5_multiple_fields(sample_dataframe):
    """Test UUID5 generation with multiple composite fields."""
    result = get_composite_uuid5(
        sample_dataframe, 
        fields=['dataset_source_id', 'dataset_title'],
        uuid_col='test_uuid'
    )
    
    assert 'test_uuid' in result.columns
    assert len(result['test_uuid'].unique()) == 3


def test_get_composite_uuid5_duplicate_detection():
    """Test that duplicate UUIDs raise an error."""
    df_with_duplicates = pd.DataFrame({
        'dataset_source_id': ['SRP123', 'SRP123'],  # Duplicate values
        'dataset_title': ['Study A', 'Study A']
    })
    
    with pytest.raises(ValueError, match="Duplicate UUID5 values found"):
        get_composite_uuid5(
            df_with_duplicates,
            fields=['dataset_source_id'],
            uuid_col='dataset_uuid'
        )


# ============================================================================
# Test PMID to SRA ID Fetching
# ============================================================================

@patch('modules.gather_sra_data.Entrez')
def test_fetch_sra_ids_success(mock_entrez, mock_entrez_response):
    """Test successful fetching of SRA IDs for a PMID."""
    mock_entrez.elink.return_value = MagicMock()
    mock_entrez.read.return_value = mock_entrez_response
    
    pmid, sra_ids, error_flag, error_msg = fetch_sra_ids('12345678')
    
    assert pmid == '12345678'
    assert sra_ids == ['SRX123456', 'SRX789012']
    assert error_flag is False
    assert error_msg == ''


@patch('modules.gather_sra_data.Entrez')
def test_fetch_sra_ids_no_results(mock_entrez):
    """Test fetching when no SRA IDs are found."""
    mock_entrez.elink.return_value = MagicMock()
    mock_entrez.read.return_value = [{'LinkSetDb': []}]
    
    pmid, sra_ids, error_flag, error_msg = fetch_sra_ids('99999999')
    
    assert pmid == '99999999'
    assert sra_ids == []
    assert error_flag is False
    assert error_msg == ''


@patch('modules.gather_sra_data.Entrez')
@patch('modules.gather_sra_data.time.sleep')
def test_fetch_sra_ids_rate_limit_retry(mock_sleep, mock_entrez):
    """Test retry logic for rate limiting errors."""
    # First call raises 429, second succeeds
    mock_entrez.elink.side_effect = [
        Exception("429 Too Many Requests"),
        MagicMock()
    ]
    mock_entrez.read.return_value = [{'LinkSetDb': [{'Link': [{'Id': 'SRX999'}]}]}]
    
    pmid, sra_ids, error_flag, error_msg = fetch_sra_ids('12345678')
    
    assert pmid == '12345678'
    assert sra_ids == ['SRX999']
    assert error_flag is False
    # Should have slept for exponential backoff
    assert mock_sleep.call_count >= 2


@patch('modules.gather_sra_data.Entrez')
def test_fetch_sra_ids_permanent_error(mock_entrez):
    """Test handling of non-rate-limiting errors."""
    mock_entrez.elink.side_effect = ValueError("Invalid PMID")
    
    pmid, sra_ids, error_flag, error_msg = fetch_sra_ids('invalid')
    
    assert pmid == 'invalid'
    assert sra_ids == []
    assert error_flag is True
    assert 'ValueError' in error_msg


@patch('modules.gather_sra_data.fetch_sra_ids')
def test_get_sra_ids_for_pubmed_ids(mock_fetch):
    """Test batch fetching of SRA IDs for multiple PMIDs."""
    # Mock responses for different PMIDs
    mock_fetch.side_effect = [
        ('12345', ['SRX123'], False, ''),
        ('67890', [], True, 'Error: timeout'),
        ('11111', ['SRX456', 'SRX789'], False, '')
    ]
    
    pmid_list = ['12345', '67890', '11111']
    results_dict, failed_dict = get_sra_ids_for_pubmed_ids(pmid_list)
    
    # Check successful results
    assert '12345' in results_dict
    assert results_dict['12345'] == ['SRX123']
    assert '11111' in results_dict
    assert results_dict['11111'] == ['SRX456', 'SRX789']
    
    # Check failed results
    assert '67890' in failed_dict
    assert failed_dict['67890'] == 'Error: timeout'


# ============================================================================
# Test Data Aggregation Functions
# ============================================================================

def test_merge_semicolon_fields():
    """Test merging and deduplication of semicolon-separated fields."""
    values = ['item1;item2', 'item2;item3', 'item1']
    result = merge_semicolon_fields(values)
    
    # Should deduplicate and sort
    assert result == 'item1; item2; item3'


def test_merge_semicolon_fields_with_empty():
    """Test merging with empty or NaN values."""
    values = ['item1;item2', '', None, 'item3']
    result = merge_semicolon_fields(values)
    
    assert result == 'item1; item2; item3'


def test_merge_semicolon_fields_all_empty():
    """Test merging when all values are empty."""
    values = ['', None, '']
    result = merge_semicolon_fields(values)
    
    assert result == ''


def test_aggregate_batch_mappings():
    """Test aggregation of multiple batch mapping DataFrames."""
    df1 = pd.DataFrame({
        'PMID': ['12345', '67890'],
        'SRA_IDs': ['SRX123', 'SRX456'],
        'SRP_ERP_IDs': ['SRP001', 'SRP002']
    })
    
    df2 = pd.DataFrame({
        'PMID': ['11111'],
        'SRA_IDs': ['SRX789'],
        'SRP_ERP_IDs': ['SRP003']
    })
    
    result = aggregate_batch_mappings([df1, df2])
    
    # Should combine all rows (no duplicates)
    assert len(result) == 3
    assert set(result['PMID']) == {'12345', '67890', '11111'}


def test_aggregate_batch_mappings_dedup_prefers_nonempty():
    """Test that dedup keeps non-empty rows over empty ones for the same PMID."""
    # First pass: PMID 12345 failed (empty SRA/SRP)
    df1 = pd.DataFrame({
        'PMID': ['12345', '67890'],
        'SRA_IDs': ['', 'SRX456'],
        'SRP_ERP_IDs': ['', 'SRP002']
    })
    
    # Second pass retry: PMID 12345 succeeded
    df2 = pd.DataFrame({
        'PMID': ['12345'],
        'SRA_IDs': ['SRX999'],
        'SRP_ERP_IDs': ['SRP001']
    })
    
    result = aggregate_batch_mappings([df1, df2])
    
    # Should deduplicate: 2 unique PMIDs, not 3 rows
    assert len(result) == 2
    assert set(result['PMID']) == {'12345', '67890'}
    
    # PMID 12345 should have the retry (non-empty) data
    row_12345 = result[result['PMID'] == '12345'].iloc[0]
    assert row_12345['SRA_IDs'] == 'SRX999'
    assert row_12345['SRP_ERP_IDs'] == 'SRP001'


def test_aggregate_batch_srp_data():
    """Test aggregation of SRP-centric batch data with deduplication."""
    df1 = pd.DataFrame({
        'SRP_ERP_ID': ['SRP001', 'SRP002'],
        'PMIDs': ['12345', '67890'],
        'PMID_Count': [1, 1]
    })
    
    df2 = pd.DataFrame({
        'SRP_ERP_ID': ['SRP001', 'SRP003'],
        'PMIDs': ['99999', '11111'],
        'PMID_Count': [1, 1]
    })
    
    result = aggregate_batch_srp_data([df1, df2])
    
    # Should group by SRP_ERP_ID and merge fields
    assert 'SRP001' in result['SRP_ERP_ID'].values
    assert 'SRP002' in result['SRP_ERP_ID'].values
    assert 'SRP003' in result['SRP_ERP_ID'].values
    
    # Check that SRP001 has merged PMIDs
    srp001_row = result[result['SRP_ERP_ID'] == 'SRP001'].iloc[0]
    assert srp001_row['PMID_Count'] == 2


# ============================================================================
# Test Batch Management Functions
# ============================================================================

def test_get_max_batch_number(sample_batch_dir):
    """Test extracting the maximum batch number from directory."""
    max_batch = get_max_batch_number(sample_batch_dir)
    assert max_batch == 2


def test_get_max_batch_number_empty_dir(tmp_path):
    """Test max batch number when directory is empty."""
    empty_dir = tmp_path / "empty_batches"
    empty_dir.mkdir()
    
    max_batch = get_max_batch_number(str(empty_dir))
    assert max_batch == 0


def test_load_processed_pmids_from_batches(sample_batch_dir):
    """Test loading PMIDs from existing batch files."""
    processed_pmids = load_processed_pmids_from_batches(sample_batch_dir)
    
    # Should include all PMIDs from both batch files
    assert '12345' in processed_pmids
    assert '67890' in processed_pmids
    assert '11111' in processed_pmids
    assert '22222' in processed_pmids
    assert len(processed_pmids) == 4


def test_load_processed_pmids_from_batches_no_dir(tmp_path):
    """Test loading PMIDs when directory doesn't exist."""
    nonexistent_dir = str(tmp_path / "nonexistent")
    
    processed_pmids = load_processed_pmids_from_batches(nonexistent_dir)
    assert processed_pmids == set()


# ============================================================================
# Test Failed PMID Reporting
# ============================================================================

@patch('modules.gather_sra_data.config')
def test_save_failed_pmids_report(mock_config, tmp_path):
    """Test saving failed PMIDs to CSV report."""
    reports_dir = tmp_path / "reports"
    mock_config.REPORTS_GATHERED_DIR = str(reports_dir)
    
    failed_dict = {
        '12345': 'Error: timeout',
        '67890': 'Error: rate limit'
    }
    
    save_failed_pmids_report(failed_dict)
    
    # Verify the file was created and contains correct data
    report_path = reports_dir / "failed_pmid_sra_searches.csv"
    assert report_path.exists(), "Report CSV was not created"
    
    saved_df = pd.read_csv(str(report_path))
    assert 'pmid' in saved_df.columns, "Missing 'pmid' column"
    assert 'failure_type' in saved_df.columns, "Missing 'failure_type' column"
    assert len(saved_df) == 2
    
    saved_dict = dict(zip(saved_df['pmid'].astype(str), saved_df['failure_type']))
    assert saved_dict['12345'] == 'Error: timeout'
    assert saved_dict['67890'] == 'Error: rate limit'


@patch('modules.gather_sra_data.config')
def test_save_failed_pmids_report_empty_dict(mock_config, tmp_path):
    """Test that saving an empty dict creates no file."""
    reports_dir = tmp_path / "reports"
    mock_config.REPORTS_GATHERED_DIR = str(reports_dir)
    
    save_failed_pmids_report({})
    
    # Should not create directory or file when dict is empty
    report_path = reports_dir / "failed_pmid_sra_searches.csv"
    assert not report_path.exists()


@patch('modules.gather_sra_data.config')
def test_load_failed_pmids_report(mock_config, tmp_path):
    """Test loading failed PMIDs from CSV report.
    
    Uses a real temporary file rather than mocking builtins.open,
    because pd.read_csv uses its own C-level file handle.
    """
    reports_dir = tmp_path / "reports"
    reports_dir.mkdir()
    mock_config.REPORTS_GATHERED_DIR = str(reports_dir)
    
    # Write a real CSV file with known content
    report_path = reports_dir / "failed_pmid_sra_searches.csv"
    report_df = pd.DataFrame({
        'pmid': ['12345', '67890'],
        'failure_type': ['timeout', 'rate_limit']
    })
    report_df.to_csv(str(report_path), index=False)
    
    failed_dict = load_failed_pmids_report()
    
    assert '12345' in failed_dict
    assert failed_dict['12345'] == 'timeout'
    assert '67890' in failed_dict
    assert failed_dict['67890'] == 'rate_limit'


@patch('modules.gather_sra_data.config')
def test_load_failed_pmids_report_no_file(mock_config, tmp_path):
    """Test loading failed PMIDs when file doesn't exist."""
    mock_config.REPORTS_GATHERED_DIR = str(tmp_path / "nonexistent_reports")
    
    failed_dict = load_failed_pmids_report()
    assert failed_dict == {}


# ============================================================================
# Test Helper Functions
# ============================================================================

def test_create_blank_metadata_dict():
    """Test creation of blank metadata dictionary."""
    result = _create_blank_metadata_dict('SRP123456')
    
    assert result['dataset_source_id'] == 'SRP123456'
    assert result['dataset_title'] == ''
    assert result['description'] == ''
    assert result['primary_disease'] == 'Not Reported'
    assert 'release_date' in result
    assert 'assay_method' in result
    assert result['study_type'] == 'Genomic sequencing'


def test_create_sra_mapping_dataframe(sample_mapping_data):
    """Test creation of PMID-to-SRA mapping DataFrame."""
    # Need to provide sra_to_srp_ids mapping as well
    sra_to_srp_ids = {
        'SRX123': ['SRP001'],
        'SRX456': ['SRP001'],
        'SRX789': ['SRP002']
    }
    
    result = create_sra_mapping_dataframe(sample_mapping_data, sra_to_srp_ids)
    
    assert 'PMID' in result.columns
    assert 'SRA_IDs' in result.columns
    assert 'SRP_ERP_IDs' in result.columns
    assert len(result) == 3  # One row per PMID
    assert '12345' in result['PMID'].values


def test_create_sra_mapping_dataframe_empty():
    """Test mapping DataFrame creation with empty input."""
    result = create_sra_mapping_dataframe({}, {})
    
    assert len(result) == 0
    assert 'PMID' in result.columns
    assert 'SRA_IDs' in result.columns
    assert 'SRP_ERP_IDs' in result.columns


def test_create_srp_centric_dataframe():
    """Test creation of SRP-centric DataFrame from mapping."""
    pmid_to_sra_ids = {
        '12345': ['SRX123', 'SRX456'],
        '67890': ['SRX789']
    }
    
    sra_to_srp_ids = {
        'SRX123': ['SRP001'],
        'SRX456': ['SRP001'],
        'SRX789': ['SRP002']
    }
    
    result = create_srp_centric_dataframe(pmid_to_sra_ids, sra_to_srp_ids)
    
    assert 'SRP_ERP_ID' in result.columns
    assert 'PMIDs' in result.columns
    assert 'PMID_Count' in result.columns
    # Should have one row per unique SRP ID
    assert len(result) == 2
    
    # Check that fields are properly semicolon-separated
    srp001_row = result[result['SRP_ERP_ID'] == 'SRP001'].iloc[0]
    assert '12345' in srp001_row['PMIDs']
    assert srp001_row['PMID_Count'] == 1


# ============================================================================
# Integration-style Tests
# ============================================================================

def test_full_uuid_generation_workflow(sample_dataframe):
    """Test complete UUID generation workflow with realistic data."""
    # This mimics the actual workflow in gather_sra_data
    result = get_composite_uuid5(
        sample_dataframe.copy(),
        fields=['dataset_source_id'],
        uuid_col='dataset_uuid'
    )
    
    # Verify all records have valid UUIDs
    assert len(result) == 3
    assert result['dataset_uuid'].notna().all()
    assert len(result['dataset_uuid'].unique()) == 3


@patch('modules.gather_sra_data.fetch_sra_ids')
def test_mapping_to_srp_workflow(mock_fetch):
    """Test the workflow from PMID list to SRP-centric data."""
    # Mock fetch responses
    mock_fetch.side_effect = [
        ('12345', ['SRX123', 'SRX456'], False, ''),
        ('67890', ['SRX789'], False, '')
    ]
    
    # Step 1: Fetch SRA IDs
    pmids = ['12345', '67890']
    results, failed = get_sra_ids_for_pubmed_ids(pmids)
    assert len(failed) == 0  # No failures expected with mocked responses
    
    # Mock SRA to SRP mapping
    sra_to_srp_ids = {
        'SRX123': ['SRP001'],
        'SRX456': ['SRP001'],
        'SRX789': ['SRP002']
    }
    
    # Step 2: Create mapping DataFrame
    mapping_df = create_sra_mapping_dataframe(results, sra_to_srp_ids)
    assert len(mapping_df) == 2  # One row per PMID
    assert 'PMID' in mapping_df.columns
    assert set(mapping_df['PMID'].values) == {'12345', '67890'}
    
    # Step 3: Create SRP-centric DataFrame
    srp_df = create_srp_centric_dataframe(results, sra_to_srp_ids)
    
    assert len(srp_df) == 2
    assert 'SRP001' in srp_df['SRP_ERP_ID'].values
    assert 'SRP002' in srp_df['SRP_ERP_ID'].values


# ============================================================================
# Tricky Edge Cases from TEST MODE (Real-world scenarios)
# ============================================================================

@patch('modules.gather_sra_data.Entrez')
def test_duplicate_pmid_handling(mock_entrez):
    """Test handling of duplicate PMIDs in input (TEST MODE case 1-2)."""
    # Mock response
    mock_entrez.elink.return_value = MagicMock()
    mock_entrez.read.return_value = [{
        'LinkSetDb': [{
            'Link': [{'Id': 'SRX123'}]
        }]
    }]
    
    # Same PMID listed twice
    pmids = ['38738472', '38738472']
    results, failed = get_sra_ids_for_pubmed_ids(pmids)
    assert len(failed) == 0  # No failures expected
    
    # Should handle duplicates gracefully
    assert '38738472' in results
    # API should only be called once due to deduplication in actual implementation
    # (or twice if no deduplication - both are valid)


@patch('modules.gather_sra_data.Entrez')
def test_no_sra_match(mock_entrez):
    """Test PMID with no SRA data (TEST MODE case 3)."""
    # Mock response with no links
    mock_entrez.elink.return_value = MagicMock()
    mock_entrez.read.return_value = [{'LinkSetDb': []}]
    
    pmid, sra_ids, error_flag, error_msg = fetch_sra_ids('10637239')
    
    assert pmid == '10637239'
    assert sra_ids == []
    assert error_flag is False  # No error, just no results
    assert error_msg == ''


@patch('modules.gather_sra_data.Entrez')
def test_erp_match(mock_entrez):
    """Test European Read Archive (ERP) ID handling (TEST MODE case 4)."""
    # Mock response with ERP IDs
    mock_entrez.elink.return_value = MagicMock()
    mock_entrez.read.return_value = [{
        'LinkSetDb': [{
            'Link': [
                {'Id': 'ERX123456'},
                {'Id': 'ERX789012'}
            ]
        }]
    }]
    
    pmid, sra_ids, error_flag, error_msg = fetch_sra_ids('26829319')
    
    assert pmid == '26829319'
    assert 'ERX123456' in sra_ids
    assert 'ERX789012' in sra_ids
    assert error_flag is False
    assert error_msg == ''


@patch('modules.gather_sra_data.Entrez')
def test_many_to_one_srp_mapping(mock_entrez):
    """Test multiple PMIDs mapping to same SRP (TEST MODE case 5-6)."""
    # Mock responses for two different PMIDs
    mock_entrez.elink.side_effect = [
        MagicMock(),
        MagicMock()
    ]
    mock_entrez.read.side_effect = [
        [{'LinkSetDb': [{'Link': [{'Id': 'SRX111'}]}]}],
        [{'LinkSetDb': [{'Link': [{'Id': 'SRX222'}]}]}]
    ]
    
    pmid1, sra_ids1, err1, _ = fetch_sra_ids('38260414')
    pmid2, sra_ids2, err2, _ = fetch_sra_ids('38802751')
    
    # Both should succeed with different SRA IDs
    assert not err1 and not err2
    assert sra_ids1 == ['SRX111']
    assert sra_ids2 == ['SRX222']
    
    # When these map to same SRP, aggregation should merge them
    pmid_to_sra = {pmid1: sra_ids1, pmid2: sra_ids2}
    sra_to_srp = {'SRX111': ['SRP999'], 'SRX222': ['SRP999']}
    
    srp_df = create_srp_centric_dataframe(pmid_to_sra, sra_to_srp)
    
    # Should have one row for the shared SRP
    assert len(srp_df) == 1
    assert srp_df.iloc[0]['SRP_ERP_ID'] == 'SRP999'
    assert srp_df.iloc[0]['PMID_Count'] == 2


@patch('modules.gather_sra_data.Entrez')
def test_bad_pmid_input(mock_entrez):
    """Test handling of invalid PMID format (TEST MODE case 7)."""
    # Mock exception for bad input
    mock_entrez.elink.side_effect = ValueError("Invalid id parameter")
    
    pmid, sra_ids, error_flag, error_msg = fetch_sra_ids('bad_input')
    
    assert pmid == 'bad_input'
    assert sra_ids == []
    assert error_flag is True
    assert 'ValueError' in error_msg


# ============================================================================
# Test SRP ID Fetching (Extended Coverage)
# ============================================================================

@patch('modules.gather_sra_data.Entrez')
def test_get_srp_ids_for_sra_id_success(mock_entrez):
    """Test fetching SRP IDs from a single SRA ID."""
    # Mock efetch to return XML with SRP accession
    mock_xml = b"""<?xml version="1.0"?>
    <EXPERIMENT_PACKAGE_SET>
        <EXPERIMENT_PACKAGE>
            <STUDY accession="SRP456">
                <DESCRIPTOR>
                    <STUDY_TITLE>Test</STUDY_TITLE>
                </DESCRIPTOR>
            </STUDY>
        </EXPERIMENT_PACKAGE>
    </EXPERIMENT_PACKAGE_SET>"""
    
    mock_handle = BytesIO(mock_xml)
    mock_entrez.efetch.return_value = mock_handle
    
    srp_ids = get_srp_ids_for_sra_id('SRX123')
    
    assert 'SRP456' in srp_ids


@patch('modules.gather_sra_data.Entrez')
def test_get_srp_ids_for_sra_id_no_results(mock_entrez):
    """Test SRA ID with no associated SRP."""
    mock_xml = b"""<?xml version="1.0"?>
    <EXPERIMENT_PACKAGE_SET>
        <EXPERIMENT_PACKAGE>
            <STUDY>
                <DESCRIPTOR>
                    <STUDY_TITLE>No accession here</STUDY_TITLE>
                </DESCRIPTOR>
            </STUDY>
        </EXPERIMENT_PACKAGE>
    </EXPERIMENT_PACKAGE_SET>"""
    
    mock_handle = BytesIO(mock_xml)
    mock_entrez.efetch.return_value = mock_handle
    
    srp_ids = get_srp_ids_for_sra_id('SRX999')
    
    assert srp_ids == []


@patch('modules.gather_sra_data.get_srp_ids_for_sra_id')
def test_get_srp_ids_for_sra_ids_batch(mock_get_srp):
    """Test batch fetching of SRP IDs for multiple SRA IDs."""
    # Mock responses
    mock_get_srp.side_effect = [
        ['SRP001'],
        ['SRP002'],
        []  # No results for third
    ]
    
    sra_ids = ['SRX123', 'SRX456', 'SRX789']
    result_dict = get_srp_ids_for_sra_ids(sra_ids)
    
    assert result_dict['SRX123'] == ['SRP001']
    assert result_dict['SRX456'] == ['SRP002']
    assert result_dict['SRX789'] == []


# ============================================================================
# Test Metadata Fetching and Parsing
# ============================================================================

@patch('modules.gather_sra_data.Entrez')
def test_fetch_sra_metadata_for_srp_success(mock_entrez, mock_xml_metadata):
    """Test successful metadata fetching and XML parsing."""
    # Mock efetch to return XML as bytes (what Entrez.efetch actually returns)
    mock_handle = BytesIO(mock_xml_metadata.encode('utf-8'))
    mock_handle.close = MagicMock()  # Mock the close method
    mock_entrez.efetch.return_value = mock_handle
    
    # Function requires both srp_id and a sample sra_id
    metadata = fetch_sra_metadata_for_srp('SRP123456', 'SRX123456')
    
    # Verify key fields were extracted
    assert metadata['dataset_source_id'] == 'SRP123456'
    assert metadata['dataset_title'] == 'Test Study Title'
    assert 'test study abstract' in metadata['description'].lower()
    assert metadata['assay_method'] == 'RNA-Seq'
    # Note: PI_name parsing requires specific XML structure not in mock
    assert metadata['release_date'] == '2023-01-15'


@patch('modules.gather_sra_data.Entrez')
def test_fetch_sra_metadata_for_srp_api_error(mock_entrez):
    """Test handling of API errors during metadata fetch."""
    mock_entrez.efetch.side_effect = Exception("API Error")
    
    metadata = fetch_sra_metadata_for_srp('SRP999', 'SRX999')
    
    # Should return blank metadata dict on error
    assert metadata['dataset_source_id'] == 'SRP999'
    assert metadata['dataset_title'] == ''
    assert metadata['description'] == ''


@patch('modules.gather_sra_data.Entrez')
def test_fetch_sra_metadata_for_srp_no_study_section(mock_entrez):
    """Test handling of valid XML that has no STUDY element."""
    xml_no_study = b"<?xml version='1.0'?><EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"
    mock_entrez.efetch.return_value = BytesIO(xml_no_study)
    
    metadata = fetch_sra_metadata_for_srp('SRP123', 'SRX123')
    
    # Should return metadata with SRP set but fields blank (no STUDY section)
    assert metadata['dataset_source_id'] == 'SRP123'
    assert metadata['dataset_title'] == ''
    assert metadata['description'] == ''
    assert metadata['dataset_source_repo'] == 'SRA'


# ============================================================================
# Test Data Enrichment
# ============================================================================

def test_enrich_srp_data_with_metadata():
    """Test enriching SRP data with metadata."""
    # Create test SRP dataframe
    srp_df = pd.DataFrame({
        'SRP_ERP_ID': ['SRP001'],
        'PMIDs': ['12345'],
        'PMID_Count': [1]
    })
    
    # Create metadata dataframe with matching column names
    metadata_df = pd.DataFrame({
        'SRP_ERP_ID': ['SRP001'],
        'dataset_source_id': ['SRP001'],
        'dataset_title': ['Test Study'],
        'description': ['Test Description'],
        'assay_method': ['RNA-Seq'],
        'PI_name': ['John Doe'],
        'release_date': ['2023-01-15'],
        'dataset_pmid': ['']
    })
    
    enriched_df = enrich_srp_data_with_metadata(srp_df, metadata_df)
    
    assert 'dataset_title' in enriched_df.columns
    assert enriched_df.iloc[0]['dataset_title'] == 'Test Study'
    assert enriched_df.iloc[0]['assay_method'] == 'RNA-Seq'
    assert enriched_df.iloc[0]['dataset_pmid'] == '12345'  # PMIDs should be merged
    # PMIDs column should be dropped (replaced by dataset_pmid)
    assert 'PMIDs' not in enriched_df.columns
    assert enriched_df.iloc[0]['PMID_Count'] == 1


@patch('modules.gather_sra_data.config')
@patch('os.path.exists')
@patch('pandas.read_csv')
def test_get_dataset_doc_from_project(mock_read_csv, mock_exists, mock_config):
    """Test deriving dataset_doc from coreproject via publication/project/program chain.
    
    The function chains: dataset_pmid → publication.pmid → publication.coreproject
    → project.project_id → project.'program.program_id' → program.program_id → program.doc
    """
    # Mock config paths
    mock_config.PUBLICATIONS_INTERMED_PATH = 'fake/pubs.csv'
    mock_config.PROJECTS_INTERMED_PATH = 'fake/projects.csv'
    mock_config.PROGRAMS_INTERMED_PATH = 'fake/programs.csv'
    
    # Mock file existence
    mock_exists.return_value = True
    
    # Mock the CSV reads with proper column names matching source code expectations
    mock_pubs = pd.DataFrame({
        'pmid': ['12345', '67890', '11111'],
        'coreproject': ['R01CA123456', 'P01CA789012', 'U01CA456789']
    })
    
    # Source code uses 'program.program_id' (dotted column name), NOT 'program_id'
    mock_projects = pd.DataFrame({
        'project_id': ['R01CA123456', 'P01CA789012', 'U01CA456789'],
        'program.program_id': ['PROG1', 'PROG2', 'PROG1']
    })
    
    # Source code uses 'doc' column, NOT 'program_doc'
    mock_programs = pd.DataFrame({
        'program_id': ['PROG1', 'PROG2'],
        'doc': ['CA', 'NCI']
    })
    
    # Setup read_csv to return different dataframes based on call order
    mock_read_csv.side_effect = [mock_pubs, mock_projects, mock_programs]
    
    # Input SRA dataframe with dataset_pmid that matches publications
    df = pd.DataFrame({
        'dataset_pmid': ['12345', '67890', '11111'],
        'dataset_source_id': ['SRP001', 'SRP002', 'SRP003']
    })
    
    result = get_dataset_doc_from_project(df)
    
    # Should have added the columns
    assert 'dataset_doc' in result.columns
    assert 'program_id' in result.columns
    assert 'funding_source' in result.columns
    assert len(result) == 3
    
    # Verify the full chain resolved correctly
    row_0 = result[result['dataset_source_id'] == 'SRP001'].iloc[0]
    assert row_0['funding_source'] == 'R01CA123456'
    assert row_0['program_id'] == 'PROG1'
    assert row_0['dataset_doc'] == 'CA'
    
    row_1 = result[result['dataset_source_id'] == 'SRP002'].iloc[0]
    assert row_1['funding_source'] == 'P01CA789012'
    assert row_1['program_id'] == 'PROG2'
    assert row_1['dataset_doc'] == 'NCI'


@patch('modules.gather_sra_data.config')
def test_get_dataset_doc_from_project_no_publication_file(mock_config, tmp_path):
    """Test dataset_doc when publication file doesn't exist."""
    mock_config.PUBLICATIONS_INTERMED_PATH = str(tmp_path / 'nonexistent_pubs.csv')
    mock_config.PROJECTS_INTERMED_PATH = str(tmp_path / 'nonexistent_projects.csv')
    mock_config.PROGRAMS_INTERMED_PATH = str(tmp_path / 'nonexistent_programs.csv')
    
    df = pd.DataFrame({
        'dataset_pmid': ['12345'],
        'dataset_source_id': ['SRP001']
    })
    
    result = get_dataset_doc_from_project(df)
    
    # Should return df_null_result with blank columns
    assert 'dataset_doc' in result.columns
    assert 'program_id' in result.columns
    assert 'funding_source' in result.columns
    assert (result['dataset_doc'] == '').all()
    assert (result['program_id'] == '').all()
    assert (result['funding_source'] == '').all()


@patch('modules.gather_sra_data.config')
@patch('os.path.exists', return_value=True)
@patch('pandas.read_csv')
def test_get_dataset_doc_from_project_missing_columns(mock_read_csv, mock_exists, mock_config):
    """Test dataset_doc when publication CSV is missing required columns."""
    mock_config.PUBLICATIONS_INTERMED_PATH = 'fake/pubs.csv'
    mock_config.PROJECTS_INTERMED_PATH = 'fake/projects.csv'
    mock_config.PROGRAMS_INTERMED_PATH = 'fake/programs.csv'
    
    # The function reads all 3 CSVs (publication, project, program) before
    # checking column requirements. Publication is missing 'coreproject'.
    mock_read_csv.side_effect = [
        pd.DataFrame({'pmid': ['12345'], 'wrong_col': ['X']}),  # pubs (missing coreproject)
        pd.DataFrame({'project_id': ['P1'], 'program.program_id': ['PR1']}),  # projects
        pd.DataFrame({'program_id': ['PR1'], 'doc': ['CA']}),  # programs
    ]
    
    df = pd.DataFrame({
        'dataset_pmid': ['12345'],
        'dataset_source_id': ['SRP001']
    })
    
    result = get_dataset_doc_from_project(df)
    
    # Should return null result due to missing 'coreproject' column
    assert 'dataset_doc' in result.columns
    assert (result['dataset_doc'] == '').all()


# ============================================================================
# Test Batch File Operations
# ============================================================================

def test_load_all_batch_files(tmp_path):
    """Test loading and combining all batch files."""
    batch_dir = tmp_path / "batches"
    batch_dir.mkdir()
    
    # Create multiple batch files
    df1 = pd.DataFrame({'col1': [1, 2], 'col2': ['a', 'b']})
    df1.to_csv(batch_dir / "batch_0001_test.csv", index=False)
    
    df2 = pd.DataFrame({'col1': [3, 4], 'col2': ['c', 'd']})
    df2.to_csv(batch_dir / "batch_0002_test.csv", index=False)
    
    # Load all files with '_test.csv' suffix
    combined = load_all_batch_files(str(batch_dir), '_test.csv')
    
    assert len(combined) == 4
    assert set(combined['col1']) == {1, 2, 3, 4}


def test_load_all_batch_files_empty_dir(tmp_path):
    """Test loading batch files from empty directory."""
    batch_dir = tmp_path / "empty_batches"
    batch_dir.mkdir()
    
    result = load_all_batch_files(str(batch_dir), '_test.csv')
    
    # Should return empty DataFrame
    assert len(result) == 0


# ============================================================================
# Test Real-world Scenarios (Integration-style)
# ============================================================================

@patch('modules.gather_sra_data.fetch_sra_ids')
def test_complete_pmid_to_enriched_datasets_workflow(mock_fetch):
    """Test the end-to-end workflow: PMIDs → SRA → SRP → mapping + metadata merge."""
    # Mock fetch_sra_ids for three PMIDs
    mock_fetch.side_effect = [
        ('12345', ['SRX100', 'SRX200'], False, ''),
        ('67890', ['SRX300'], False, ''),
        ('11111', [], False, ''),  # No SRA match
    ]
    
    # Step 1: Fetch SRA IDs
    pmids = ['12345', '67890', '11111']
    pmid_to_sra, failed = get_sra_ids_for_pubmed_ids(pmids)
    assert len(failed) == 0
    assert pmid_to_sra['12345'] == ['SRX100', 'SRX200']
    assert pmid_to_sra['11111'] == []
    
    # Step 2: Simulate SRA → SRP mapping
    sra_to_srp = {
        'SRX100': ['SRP001'],
        'SRX200': ['SRP001'],  # Same SRP as SRX100
        'SRX300': ['SRP002'],
    }
    
    # Step 3a: Create mapping DataFrame
    mapping_df = create_sra_mapping_dataframe(pmid_to_sra, sra_to_srp)
    assert len(mapping_df) == 3  # One row per PMID
    row_12345 = mapping_df[mapping_df['PMID'] == '12345'].iloc[0]
    assert 'SRP001' in row_12345['SRP_ERP_IDs']
    
    # Step 3b: Create SRP-centric DataFrame
    srp_df = create_srp_centric_dataframe(pmid_to_sra, sra_to_srp)
    assert len(srp_df) == 2  # SRP001 and SRP002
    
    # Step 4: Create SRP→SRA mapping for metadata
    srp_to_sra = create_srp_to_sra_mapping(pmid_to_sra, sra_to_srp)
    assert 'SRP001' in srp_to_sra
    assert 'SRP002' in srp_to_sra
    
    # Step 5: Merge with metadata
    metadata_df = pd.DataFrame({
        'SRP_ERP_ID': ['SRP001', 'SRP002'],
        'dataset_source_id': ['SRP001', 'SRP002'],
        'dataset_title': ['Study A', 'Study B'],
        'description': ['Desc A', 'Desc B'],
        'dataset_pmid': ['', ''],
    })
    
    enriched = enrich_srp_data_with_metadata(srp_df, metadata_df)
    assert len(enriched) == 2
    assert 'dataset_title' in enriched.columns
    # PMIDs from srp_df should merge into dataset_pmid
    srp001_row = enriched[enriched['SRP_ERP_ID'] == 'SRP001'].iloc[0]
    assert '12345' in srp001_row['dataset_pmid']


@patch('modules.gather_sra_data.fetch_sra_ids')
@patch('modules.gather_sra_data.get_srp_ids_for_sra_ids')
def test_workflow_with_mixed_results(mock_get_srp, mock_fetch_sra):
    """Test workflow with mix of successes, failures, and empty results."""
    # Mock diverse responses
    mock_fetch_sra.side_effect = [
        ('12345', ['SRX123'], False, ''),      # Success
        ('67890', [], False, ''),               # No results
        ('99999', [], True, 'API Error'),       # Error
    ]
    
    mock_get_srp.return_value = {
        'SRX123': ['SRP001']
    }
    
    pmids = ['12345', '67890', '99999']
    results, failed = get_sra_ids_for_pubmed_ids(pmids)
    
    # Check successful result
    assert '12345' in results
    assert results['12345'] == ['SRX123']
    
    # Check empty result (not a failure)
    assert '67890' in results
    assert results['67890'] == []
    
    # Check actual failure
    assert '99999' in failed
    assert 'API Error' in failed['99999']


# ============================================================================
# Test create_srp_to_sra_mapping
# ============================================================================

def test_create_srp_to_sra_mapping_basic():
    """Test basic SRP to sample SRA mapping creation."""
    pmid_to_sra_ids = {
        '12345': ['SRX100', 'SRX200'],
        '67890': ['SRX300'],
    }
    sra_to_srp_ids = {
        'SRX100': ['SRP001'],
        'SRX200': ['SRP001', 'SRP002'],
        'SRX300': ['SRP003'],
    }
    
    result = create_srp_to_sra_mapping(pmid_to_sra_ids, sra_to_srp_ids)
    
    # Should map each SRP to ONE sample SRA (the first encountered)
    assert 'SRP001' in result
    assert 'SRP002' in result
    assert 'SRP003' in result
    # SRP001 should get the first SRA that maps to it
    assert result['SRP001'] == 'SRX100'
    assert result['SRP003'] == 'SRX300'


def test_create_srp_to_sra_mapping_empty():
    """Test SRP to SRA mapping with empty inputs."""
    result = create_srp_to_sra_mapping({}, {})
    assert result == {}


def test_create_srp_to_sra_mapping_no_srp_matches():
    """Test when SRA IDs have no SRP matches."""
    pmid_to_sra_ids = {'12345': ['SRX100']}
    sra_to_srp_ids = {}  # No SRA→SRP mappings
    
    result = create_srp_to_sra_mapping(pmid_to_sra_ids, sra_to_srp_ids)
    assert result == {}


# ============================================================================
# Test gather_srp_metadata
# ============================================================================

@patch('modules.gather_sra_data.fetch_sra_metadata_for_srp')
def test_gather_srp_metadata_basic(mock_fetch):
    """Test gathering metadata for a list of SRP IDs."""
    # Mock the per-SRP fetch
    mock_fetch.side_effect = [
        {
            'dataset_source_id': 'SRP001',
            'dataset_title': 'Study A',
            'description': 'Desc A',
            'type': 'dataset',
            'dataset_source_repo': 'SRA',
        },
        {
            'dataset_source_id': 'SRP002',
            'dataset_title': 'Study B',
            'description': 'Desc B',
            'type': 'dataset',
            'dataset_source_repo': 'SRA',
        },
    ]
    
    srp_ids = ['SRP001', 'SRP002']
    srp_to_sra = {'SRP001': 'SRX100', 'SRP002': 'SRX200'}
    
    result = gather_srp_metadata(srp_ids, srp_to_sra, max_workers=1)
    
    assert 'SRP_ERP_ID' in result.columns
    assert len(result) == 2
    assert set(result['SRP_ERP_ID']) == {'SRP001', 'SRP002'}


@patch('modules.gather_sra_data.fetch_sra_metadata_for_srp')
def test_gather_srp_metadata_missing_sra_sample(mock_fetch):
    """Test gather_srp_metadata when SRP has no SRA sample in the mapping."""
    # SRP002 has no SRA sample → function should call _create_blank_metadata_dict
    mock_fetch.return_value = {
        'dataset_source_id': 'SRP001',
        'dataset_title': 'Study A',
        'type': 'dataset',
        'dataset_source_repo': 'SRA',
    }
    
    srp_ids = ['SRP001', 'SRP002']
    srp_to_sra = {'SRP001': 'SRX100'}  # SRP002 not mapped
    
    result = gather_srp_metadata(srp_ids, srp_to_sra, max_workers=1)
    
    assert len(result) == 2
    # SRP002 should have blank metadata
    srp002_row = result[result['dataset_source_id'] == 'SRP002']
    assert len(srp002_row) == 1
    assert srp002_row.iloc[0]['dataset_title'] == ''


# ============================================================================
# Test JSON Batch Persistence (sra_to_srp_ids)
# ============================================================================

def test_json_batch_roundtrip(tmp_path):
    """Test that sra_to_srp_ids persisted as JSON can be loaded accurately."""
    batch_dir = tmp_path / "sra_batches"
    batch_dir.mkdir()
    
    # Simulate what gather_sra_data does: save per-batch JSON
    sra_to_srp = {
        'SRX100': ['SRP001', 'SRP002'],
        'SRX200': ['SRP001'],
        'SRX300': [],
    }
    json_path = batch_dir / "batch_0001_sra_to_srp.json"
    with open(str(json_path), 'w', encoding='utf-8') as f:
        json.dump(sra_to_srp, f)
    
    # Load it back
    with open(str(json_path), 'r', encoding='utf-8') as f:
        loaded = json.load(f)
    
    assert loaded == sra_to_srp
    assert loaded['SRX100'] == ['SRP001', 'SRP002']
    assert loaded['SRX300'] == []


def test_json_batch_preserves_per_sra_granularity(tmp_path):
    """Test that JSON preserves per-SRA→SRP mapping that CSV aggregation loses.
    
    This validates the fix for the reconstruction bug: the CSV aggregates
    SRP IDs per PMID, losing per-SRA granularity. The JSON retains it.
    """
    batch_dir = tmp_path / "sra_batches"
    batch_dir.mkdir()
    
    # Two SRAs from same PMID map to different SRPs
    pmid_to_sra = {'12345': ['SRX100', 'SRX200']}
    sra_to_srp = {
        'SRX100': ['SRP001'],
        'SRX200': ['SRP002'],
    }
    
    # Save mapping CSV (aggregated per PMID - loses granularity)
    mapping_df = create_sra_mapping_dataframe(pmid_to_sra, sra_to_srp)
    mapping_path = batch_dir / "batch_0001_mapping.csv"
    mapping_df.to_csv(str(mapping_path), index=False)
    
    # Save JSON (preserves per-SRA granularity)
    json_path = batch_dir / "batch_0001_sra_to_srp.json"
    with open(str(json_path), 'w', encoding='utf-8') as f:
        json.dump(sra_to_srp, f)
    
    # Load CSV - shows both SRPs aggregated for PMID 12345
    loaded_csv = pd.read_csv(str(mapping_path))
    srp_cell = loaded_csv.iloc[0]['SRP_ERP_IDs']
    # CSV has both SRPs in one cell, can't tell which SRA goes with which SRP
    assert 'SRP001' in srp_cell
    assert 'SRP002' in srp_cell
    
    # Load JSON - preserves exact SRA→SRP mapping
    with open(str(json_path), 'r', encoding='utf-8') as f:
        loaded_json = json.load(f)
    assert loaded_json['SRX100'] == ['SRP001']
    assert loaded_json['SRX200'] == ['SRP002']


# ============================================================================
# Test get_max_batch_number Edge Cases
# ============================================================================

def test_get_max_batch_number_nonexistent_dir():
    """Test max batch number when directory doesn't exist at all."""
    max_batch = get_max_batch_number('/nonexistent/path/that/cannot/exist')
    assert max_batch == 0


def test_get_max_batch_number_no_matching_files(tmp_path):
    """Test max batch number when directory has files but none match pattern."""
    batch_dir = tmp_path / "batches"
    batch_dir.mkdir()
    
    # Create files that don't match the batch_*_mapping.csv pattern
    (batch_dir / "random_file.csv").write_text("data")
    (batch_dir / "batch_notes.txt").write_text("notes")
    
    max_batch = get_max_batch_number(str(batch_dir))
    assert max_batch == 0


# ============================================================================
# Test enrich_srp_data_with_metadata edge cases
# ============================================================================

def test_enrich_srp_data_merges_metadata_and_mapping_pmids():
    """Test that PMIDs from both mapping and metadata are merged/deduped."""
    srp_df = pd.DataFrame({
        'SRP_ERP_ID': ['SRP001'],
        'PMIDs': ['12345; 67890'],
        'PMID_Count': [2]
    })
    
    metadata_df = pd.DataFrame({
        'SRP_ERP_ID': ['SRP001'],
        'dataset_source_id': ['SRP001'],
        'dataset_title': ['Study'],
        'dataset_pmid': ['67890; 99999'],  # 67890 overlaps, 99999 is new
    })
    
    enriched = enrich_srp_data_with_metadata(srp_df, metadata_df)
    
    # dataset_pmid should contain merged, deduped PMIDs
    pmid_str = enriched.iloc[0]['dataset_pmid']
    pmids = [p.strip() for p in pmid_str.split(';')]
    assert '12345' in pmids
    assert '67890' in pmids
    assert '99999' in pmids
    assert len(pmids) == 3  # No duplicates


# ============================================================================
# Test aggregate_batch_mappings edge cases
# ============================================================================

def test_aggregate_batch_mappings_empty_input():
    """Test aggregation with empty DataFrame list."""
    empty_df = pd.DataFrame(columns=['PMID', 'SRA_IDs', 'SRP_ERP_IDs'])
    result = aggregate_batch_mappings([empty_df])
    assert len(result) == 0


def test_aggregate_batch_mappings_no_pmid_column():
    """Test aggregation when DataFrame has no PMID column (early return)."""
    df = pd.DataFrame({'other_col': ['a', 'b']})
    result = aggregate_batch_mappings([df])
    # Should return combined df without dedup processing
    assert len(result) == 2


if __name__ == "__main__":

    pytest.main([__file__, "-v"])
