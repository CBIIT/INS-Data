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
    detect_invalid_nofos,
    detect_invalid_awards,
    create_program_id,
    validate_nofos,
    validate_awards,
    report_invalid_nofos,
    report_invalid_awards,
    validate_and_rename_columns,
    replace_blank_values,
    force_replace_comma_separation,
    remove_extra_default_list_values,
    check_for_duplicate_names,
    generate_program_id_column,
    apply_value_fix,
    prompt_to_continue
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


def test_clean_nofo_and_award_cols_not_found():
    """Test cleaning NOFO and Award columns error handling when cols not found."""
    missing_cols_df = pd.DataFrame({'Column1': ['value1','value2']})

    with pytest.raises(ValueError,
                    match=f"Columns 'nofo' and 'award' must be present in the DataFrame"):
        clean_nofo_and_award_cols(missing_cols_df)


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


@pytest.fixture
def sample_bad_nofo_award_dataframe():
    return pd.DataFrame({
        'program_name': ['Program1', 'Program2'],
        'nofo': ['bad_nofo1', 'bad_nofo2; bad_nofo3'],
        'award': ['bad_award1', 'bad_award2; bad_award3'],
    })


def test_detect_invalid_nofos(sample_bad_nofo_award_dataframe):
    """Test generation of invalid nofo df from bad progam df."""
    invalid_nofos_df = detect_invalid_nofos(sample_bad_nofo_award_dataframe)
    
    assert len(invalid_nofos_df) > 0
    assert isinstance(invalid_nofos_df, pd.DataFrame)
    assert 'program_name' in invalid_nofos_df.columns
    assert 'invalid_nofo' in invalid_nofos_df.columns


def test_detect_invalid_nofos_none_found(sample_dataframe):
    """Test generation of invalid nofo df from bad progam df."""
    invalid_nofos_df = detect_invalid_nofos(sample_dataframe)
    
    assert len(invalid_nofos_df) == 0


def test_detect_invalid_awards(sample_bad_nofo_award_dataframe):
    """Test generation of invalid award df from bad progam df."""
    invalid_awards_df = detect_invalid_awards(sample_bad_nofo_award_dataframe)
    
    assert len(invalid_awards_df) > 0
    assert isinstance(invalid_awards_df, pd.DataFrame)
    assert 'program_name' in invalid_awards_df.columns
    assert 'invalid_award' in invalid_awards_df.columns


def test_detect_invalid_awards_none_found(sample_dataframe):
    """Test generation of invalid award df from bad progam df."""
    invalid_awards_df = detect_invalid_awards(sample_dataframe)
    
    assert len(invalid_awards_df) == 0


def test_create_program_id():
    """Test program_id generation from various inputs."""
    assert create_program_id('Test Program') == 'test_program'
    assert create_program_id('Test-Program') == 'test_program'
    assert create_program_id('Test Program!') == 'test_program'
    assert create_program_id('UPPER Case') == 'upper_case'
    assert create_program_id('Non-Alpha & Num*^()Char') == 'non_alpha_num_char'


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


def test_report_invalid_nofos(tmp_path, sample_bad_nofo_award_dataframe):
    """Test export of invalid nofo list"""
    report_path = tmp_path / 'invalid_nofos_report.csv'
    
    invalid_nofos_df = detect_invalid_nofos(sample_bad_nofo_award_dataframe)
    report_invalid_nofos(invalid_nofos_df, str(report_path), printout=True)
    
    assert os.path.exists(str(report_path))
    reported_df = pd.read_csv(str(report_path))
    assert not reported_df.empty


def test_report_invalid_nofos_none_found(tmp_path, sample_dataframe):
    """Test export of invalid nofo list when none present"""
    report_path = tmp_path / 'invalid_nofos_report.csv'
    
    invalid_nofos_df = detect_invalid_nofos(sample_dataframe)
    report_invalid_nofos(invalid_nofos_df, str(report_path), printout=True)
    
    assert invalid_nofos_df.empty
    assert not os.path.exists(str(report_path))


def test_report_invalid_awards(tmp_path, sample_bad_nofo_award_dataframe):
    """Test export of invalid award list"""
    report_path = tmp_path / 'invalid_awards_report.csv'
    
    invalid_awards_df = detect_invalid_awards(sample_bad_nofo_award_dataframe)
    report_invalid_awards(invalid_awards_df, str(report_path), printout=True)
    
    assert os.path.exists(str(report_path))
    reported_df = pd.read_csv(str(report_path))
    assert not reported_df.empty


def test_report_invalid_awards_none_found(tmp_path, sample_dataframe):
    """Test export of invalid award list when none present"""
    report_path = tmp_path / 'invalid_awards_report.csv'
    
    invalid_awards_df = detect_invalid_awards(sample_dataframe)
    report_invalid_awards(invalid_awards_df, str(report_path), printout=True)
    
    assert invalid_awards_df.empty
    assert not os.path.exists(str(report_path))


def test_validate_and_rename_columns():
    """Test the column validation and renaming with good input."""
    df = pd.DataFrame({
        'Old Column 1': [1, 2, 3],
        'Old Column 2': ['A', 'B', 'C'],
        'Old Column 3\r\n(with newline)': [4, 5, 6],
        'login_id': ['user1', 'user2', 'user3']
    })
    
    col_dict = {
        'Old Column 1': 'new_col1',
        'Old Column 2': 'new_col2',
        'Old Column 3\n(with newline)': 'new_col3',
        'login_id': 'login_id'
    }
    
    renamed_df = validate_and_rename_columns(df, col_dict)
    
    assert list(renamed_df.columns) == ['new_col1', 'new_col2', 'new_col3']


def test_validate_and_rename_columns_mismatch():
    """Test the column validation and renaming with bad input for error handling."""
    df = pd.DataFrame({
        'Wrong Column 1': [1, 2, 3],
        'Wrong Column 2': ['A', 'B', 'C']
    })
    
    col_dict = {
        'Correct Column 1': 'new_col1',
        'Correct Column 2': 'new_col2'
    }
    
    with pytest.raises(ValueError, match='Column names do not match expected'):
        validate_and_rename_columns(df, col_dict)


def test_validate_and_rename_columns_last_col():
    """Test the column validation and renaming last column error handling"""
    df = pd.DataFrame({
        'Old Column 1': [1, 2, 3],
        'Old Column 2': ['A', 'B', 'C'],
    })

    col_dict = {
        'Old Column 1': 'new_col1',
        'Old Column 2': 'new_col2',
    }

    with pytest.raises(ValueError, match='Unexpected final column'):
        validate_and_rename_columns(df, col_dict)


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
        'program_name': ['Program1', 'Program, comma'],
        'focus_area': ['Area1,Area2', 'Area3, Area4'],
        'cancer_type': ['Type1,Type2', 'Type3, Type4']
    })
    
    df_replaced = force_replace_comma_separation(df, ['focus_area', 
                                                      'cancer_type',
                                                      'invalid_col'])

    assert df_replaced['program_name'].tolist() == ['Program1', 'Program, comma']
    assert df_replaced['focus_area'].tolist() == ['Area1;Area2', 'Area3; Area4']
    assert df_replaced['cancer_type'].tolist() == ['Type1;Type2', 'Type3; Type4']


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


@pytest.fixture
def sample_df_with_duplicates():
    """Fixture to create a DataFrame with duplicate names and acronyms."""
    return pd.DataFrame({
        'program_name': ['Duplicate Program', 'Duplicate Program', 'Unique Program'],
        'program_acronym': ['DUP1', 'DUP1', 'UNIQ']
    })


@pytest.fixture
def sample_df_without_duplicates():
    """Fixture to create a DataFrame without duplicates."""
    return pd.DataFrame({
        'program_name': ['Unique Program 1', 'Unique Program 2'],
        'program_acronym': ['UNIQ1', 'UNIQ2']
    })


def test_check_for_duplicate_names_with_duplicates_continue(sample_df_with_duplicates, 
                                                            monkeypatch):
    """Test check_for_duplicate_names with duplicate names and user continue."""
    # Simulate user input to continue
    monkeypatch.setattr('builtins.input', lambda _: 'y')
    
    result = check_for_duplicate_names(sample_df_with_duplicates)
    assert result is True


def test_check_for_duplicate_names_with_duplicates_abort(sample_df_with_duplicates,
                                                         monkeypatch):
    """Test check_for_duplicate_names with duplicate names and user abort."""
    # Simulate user input to continue
    monkeypatch.setattr('builtins.input', lambda _: 'n')
    
    result = check_for_duplicate_names(sample_df_with_duplicates)
    assert result is False


def test_check_for_duplicate_names_without_duplicates(sample_df_without_duplicates):
    """Test check_for_duplicate_names without duplicates."""
    result = check_for_duplicate_names(sample_df_without_duplicates)
    assert result is True


def test_generate_program_id_column():
    """Test generating program_id column with various scenarios."""
    df = pd.DataFrame({
        'program_acronym': ['PROG1', '', 'NAN'],
        'program_name': ['First Program', 'Second Program', 'Third Program']
    })
    
    df_with_ids = generate_program_id_column(df)
    
    expected_ids = ['prog1', 'second_program', 'third_program']
    assert df_with_ids['program_id'].tolist() == expected_ids


def test_apply_value_fix():
    """Test applying fixes to a DataFrame."""
    key_programs_df = pd.DataFrame({
        'program_name': ['Prog A', 'Prog B'],
        'nofo': ['BAD-123', 'GOOD-456']
    })
    
    reviewed_df = pd.DataFrame({
        'program_name': ['Prog A'],
        'invalid_nofo': ['BAD-123'],
        'suggested_fix': ['GOOD-123']
    })
    
    updated_df = apply_value_fix(key_programs_df, reviewed_df, 'nofo')
    
    assert updated_df.loc[0, 'nofo'] == 'GOOD-123'
    assert updated_df.loc[1, 'nofo'] == 'GOOD-456'


def test_prompt_to_continue_empty_df(monkeypatch):
    """Test prompt_to_continue with an empty DataFrame."""
    empty_df = pd.DataFrame()
    
    result = prompt_to_continue(empty_df)
    assert result is True


@pytest.fixture
def sample_invalid_df_for_prompt():
    """Fixture to create a DataFrame with invalid values for the continue prompt."""
    return pd.DataFrame({
        'program_name': ['Program1', 'Program2'],
        'invalid_nofo': ['invalid_val1', 'invalid_val2'],
        })


def test_prompt_to_continue_non_empty_df(sample_invalid_df_for_prompt, monkeypatch):
    """Test prompt_to_continue with a non-empty DataFrame."""
    # Simulate user input for continuing
    monkeypatch.setattr('builtins.input', lambda _: 'Y')
    
    result = prompt_to_continue(sample_invalid_df_for_prompt)
    assert result is True


def test_prompt_to_continue_user_abort(sample_invalid_df_for_prompt, monkeypatch):
    """Test prompt_to_continue when user chooses to abort."""
    # Simulate user input to abort
    monkeypatch.setattr('builtins.input', lambda _: 'N')
    
    result = prompt_to_continue(sample_invalid_df_for_prompt)
    assert result is False