"""
test_gather_program_data.py
2024-12-11 ZD

Pytest test suite for the `gather_program_data.py` module.
"""

import os
import sys
import pandas as pd
import pytest

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.gather_program_data import (
    find_header_location,
    drop_obsolete_columns,
    clean_nofo_and_award_cols,
    clean_and_split_nofo_award_strings,
    is_valid_nofo,
    is_valid_award,
    create_program_id,
    validate_nofos,
    validate_awards,
    replace_blank_values,
    force_replace_comma_separation,
    remove_extra_default_list_values
)


@pytest.fixture
def test_qualtrics_csv():
    """Fixture to provide path to test Qualtrics CSV."""
    test_file_path = os.path.join(os.path.dirname(__file__), 
                                  "test_files", "test_qualtrics_output.csv")
    return test_file_path


@pytest.fixture
def sample_dataframe():
    """Fixture to create a sample DataFrame for testing."""
    return pd.DataFrame({
        'program_name': ['Test Program 1', 'Test Program 2'],
        'nofo': ['PAR-22-123', 'RFA-CA-22-456'],
        'award': ['1R01CA123456', '2P01CA789012'],
        'obsolete_column1': ['value1', 'value2'],
        'obsolete_column2': ['old1', 'old2']
    })


def test_find_header_location_success(test_qualtrics_csv):
    """Test finding header location in CSV with valid key."""
    key_value = "Name of Key Program"
    row, col = find_header_location(test_qualtrics_csv, key_value)
    
    raw_df = pd.read_csv(test_qualtrics_csv, header=None)
    test_key_value = raw_df.iloc[row, col]
    assert test_key_value == key_value


def test_find_header_location_fail(test_qualtrics_csv):
    """Test finding header location with invalid key raises ValueError."""
    key_value = "invalid_value_for_testing"
    with pytest.raises(ValueError,
                       match=f"Key value '{key_value}' not found in file."):
        find_header_location(test_qualtrics_csv, key_value)


def test_drop_obsolete_columns(sample_dataframe):
    """Test dropping columns containing 'obsolete' string."""
    df_cleaned = drop_obsolete_columns(sample_dataframe)
    
    assert 'obsolete_column1' not in df_cleaned.columns
    assert 'obsolete_column2' not in df_cleaned.columns
    assert 'program_name' in df_cleaned.columns


def test_clean_nofo_and_award_cols(sample_dataframe):
    """Test cleaning NOFO and Award columns."""
    df_cleaned = clean_nofo_and_award_cols(sample_dataframe)
    
    assert df_cleaned['nofo'].tolist() == ['PAR-22-123', 'RFA-CA-22-456']
    assert df_cleaned['award'].tolist() == ['1R01CA123456', '2P01CA789012']


def test_clean_and_split_nofo_award_strings():
    """Test parsing string list into actual lists."""
    assert clean_and_split_nofo_award_strings('PAR-22-123; RFA-CA-22-456') == \
           ['PAR-22-123', 'RFA-CA-22-456']
    
    assert clean_and_split_nofo_award_strings('') == []
    assert clean_and_split_nofo_award_strings('nan') == []


@pytest.mark.parametrize("nofo,expected", [
    ('RFA-CA-22-123', True),
    ('PAR-22-123', True),
    ('OTA-22-123', True),
    ('CA-22-123', True),
    ('invalid-nofo', False),
    ('RFA22-123', False),
    ('123', False)
])

def test_is_valid_nofo(nofo, expected):
    """Test NOFO validation with various input formats."""
    assert is_valid_nofo(nofo) == expected


@pytest.mark.parametrize("award,expected", [
    ('1R01CA123456', True),
    ('1R01CA123456-01', True),
    ('1R01CA123456-01A1', True),
    ('1R01CA123456-01S1', True),
    ('123456ABCDEFGHIJKLMNOP', True),
    ('invalid-award', False),
    ('1234', False)
])

def test_is_valid_award(award, expected):
    """Test Award validation with various input formats."""
    assert is_valid_award(award) == expected


def test_create_program_id():
    """Test program_id generation from various inputs."""
    assert create_program_id('Test Program') == 'test_program'
    assert create_program_id('Test-Program') == 'test_program'
    assert create_program_id('Test Program!') == 'test_program'
    assert create_program_id('UPPER Case') == 'upper_case'


def test_validate_nofos():
    """Test NOFO validation with good and bad inputs."""
    good_nofos = ['PAR-22-123', 'RFA-CA-22-456']
    bad_nofos = ['invalid-nofo', 'bad-format']
    
    good_result = validate_nofos('Test Program', good_nofos)
    assert good_result.empty
    
    bad_result = validate_nofos('Test Program', bad_nofos)
    assert not bad_result.empty
    assert len(bad_result) == 2
    assert all(bad_result['program_name'] == 'Test Program')


def test_validate_awards():
    """Test Award validation with good and bad inputs."""
    good_awards = ['1R01CA123456', '2P01CA789012']
    bad_awards = ['invalid-award', 'short']
    
    good_result = validate_awards('Test Program', good_awards)
    assert good_result.empty
    
    bad_result = validate_awards('Test Program', bad_awards)
    assert not bad_result.empty
    assert len(bad_result) == 2
    assert all(bad_result['program_name'] == 'Test Program')


def test_replace_blank_values():
    """Test replacing blank values in DataFrame."""
    df = pd.DataFrame({
        'column1': [None, 'value', ''],
        'column2': [None, 'data', ''],
        'column3': [None, 'white space', '']
        })
    
    blank_dict = {
        'column1': 'Default1',
        'column2': 'Default2',
        'column3': 'Default 3'
    }
    
    df_replaced = replace_blank_values(df, blank_dict)
    
    assert df_replaced['column1'].tolist() == ['Default1', 'value', 'Default1']
    assert df_replaced['column2'].tolist() == ['Default2', 'data', 'Default2']
    assert df_replaced['column3'].tolist() == ['Default 3', 'white space', 'Default 3']


def test_force_replace_comma_separation():
    """Test replacing commas with semicolons in specified columns."""
    df = pd.DataFrame({
        'focus_area': ['Area1,Area2', 'Area3,Area4'],
        'cancer_type': ['Type1,Type2', 'Type3,Type4']
    })
    
    df_replaced = force_replace_comma_separation(df, ['focus_area', 'cancer_type'])
    
    assert df_replaced['focus_area'].tolist() == ['Area1;Area2', 'Area3;Area4']
    assert df_replaced['cancer_type'].tolist() == ['Type1;Type2', 'Type3;Type4']


def test_remove_extra_default_list_values():
    """Test removing default values when other values are present."""
    df = pd.DataFrame({
        'focus_area': ['Multiple Areas;Specific Area', 
                       'Specific Area', 
                       'Multiple Areas'],
        'program_acronym': ['PROG1', 'PROG2', 'PROG3']
    })
    
    blank_dict = {
        'focus_area': 'Multiple Areas'
    }
    
    df_cleaned = remove_extra_default_list_values(df, blank_dict)
    
    expected_focus_areas = [
        'Specific Area', 
        'Specific Area', 
        'Multiple Areas'
    ]
    
    assert df_cleaned['focus_area'].tolist() == expected_focus_areas