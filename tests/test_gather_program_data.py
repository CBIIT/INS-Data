"""
test_gather_program_data.py
2024-12-11 ZD

Pytest test suite for the `gather_program_data.py` module.
"""

import os
import pandas as pd
import pytest
from unittest.mock import patch

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
    prompt_to_continue,
    load_and_clean_programs,
    gather_program_data,
)


# ============================================================
# Shared Fixtures
# ============================================================

@pytest.fixture
def test_qualtrics_csv():
    """Path to test Qualtrics CSV with real data."""
    return os.path.join(os.path.dirname(__file__),
                        "test_files", "test_qualtrics_output.csv")


@pytest.fixture
def sample_dataframe():
    """Sample DataFrame with valid NOFOs, Awards, and obsolete columns."""
    return pd.DataFrame({
        'program_name': ['Test Program 1', 'Test Program 2'],
        'nofo': ['PAR-22-123', 'RFA-CA-22-456'],
        'award': ['1R01CA123456', '2P01CA789012'],
        'obsolete_column1': ['value1', 'value2'],
        'obsolete_column2': ['old1', 'old2'],
    })


@pytest.fixture
def valid_program_csv():
    """Path to a simplified valid programs input CSV."""
    return os.path.join(os.path.dirname(__file__),
                        "test_files", "programs_valid_sample.csv")


@pytest.fixture
def sample_column_dict():
    """Standard column mapping for the simplified test CSV."""
    return {
        'Program Name': 'program_name',
        'Program Acronym': 'program_acronym',
        'NOFO': 'nofo',
        'Award': 'award',
        'Focus Area': 'focus_area',
        'Cancer Type': 'cancer_type',
        'login_id': 'login_id',
    }


# ============================================================
# CSV Parsing
# ============================================================

class TestFindHeaderLocation:

    def test_success(self, test_qualtrics_csv):
        """Finds the correct cell for a known header key."""
        key_value = "Name of Key Program"
        row, col = find_header_location(test_qualtrics_csv, key_value)

        raw_df = pd.read_csv(test_qualtrics_csv, header=None)
        assert raw_df.iloc[row, col] == key_value

    def test_key_not_found(self, test_qualtrics_csv):
        """Raises ValueError when the key doesn't exist in the file."""
        key_value = "invalid_value_for_testing"
        with pytest.raises(ValueError,
                           match=f"Key value '{key_value}' not found in file."):
            find_header_location(test_qualtrics_csv, key_value)


class TestDropObsoleteColumns:

    def test_drops_obsolete_keeps_others(self, sample_dataframe):
        df_cleaned = drop_obsolete_columns(sample_dataframe)

        assert 'obsolete_column1' not in df_cleaned.columns
        assert 'obsolete_column2' not in df_cleaned.columns
        assert 'program_name' in df_cleaned.columns


# ============================================================
# NOFO / Award Cleaning and Validation
# ============================================================

class TestCleanNofoAndAwardCols:

    def test_cleans_existing_columns(self, sample_dataframe):
        df_cleaned = clean_nofo_and_award_cols(sample_dataframe)

        assert df_cleaned['nofo'].tolist() == ['PAR-22-123', 'RFA-CA-22-456']
        assert df_cleaned['award'].tolist() == ['1R01CA123456', '2P01CA789012']

    def test_missing_columns_raises(self):
        df = pd.DataFrame({'Column1': ['value1', 'value2']})
        with pytest.raises(ValueError,
                           match="Columns 'nofo' and 'award' must be present"):
            clean_nofo_and_award_cols(df)


class TestCleanAndSplitNofoAwardStrings:

    def test_splits_semicolons(self):
        assert clean_and_split_nofo_award_strings('PAR-22-123; RFA-CA-22-456') == \
               ['PAR-22-123', 'RFA-CA-22-456']

    def test_empty_and_nan(self):
        assert clean_and_split_nofo_award_strings('') == []
        assert clean_and_split_nofo_award_strings('nan') == []


class TestIsValidNofo:

    @pytest.mark.parametrize("nofo,expected", [
        ('RFA-CA-22-123', True),
        ('PAR-22-123', True),
        ('OTA-22-123', True),
        ('CA-22-123', True),
        ('invalid-nofo', False),
        ('RFA22-123', False),
        ('123', False),
    ])
    def test_validation(self, nofo, expected):
        assert is_valid_nofo(nofo) == expected


class TestIsValidAward:

    @pytest.mark.parametrize("award,expected", [
        ('1R01CA123456', True),
        ('1R01CA123456-01', True),
        ('1R01CA123456-01A1', True),
        ('1R01CA123456-01S1', True),
        ('123456ABCDEFGHIJKLMNOP', True),
        ('invalid-award', False),
        ('1234', False),
    ])
    def test_validation(self, award, expected):
        assert is_valid_award(award) == expected


class TestValidateNofos:

    def test_good_input_returns_empty(self):
        result = validate_nofos('Test Program', ['PAR-22-123', 'RFA-CA-22-456'])
        assert result.empty

    def test_bad_input_returns_invalid_rows(self):
        result = validate_nofos('Test Program', ['invalid-nofo', 'bad-format'])

        assert len(result) == 2
        assert all(result['program_name'] == 'Test Program')


class TestValidateAwards:

    def test_good_input_returns_empty(self):
        result = validate_awards('Test Program', ['1R01CA123456', '2P01CA789012'])
        assert result.empty

    def test_bad_input_returns_invalid_rows(self):
        result = validate_awards('Test Program', ['invalid-award', 'short'])

        assert len(result) == 2
        assert all(result['program_name'] == 'Test Program')


# ============================================================
# Invalid NOFO / Award Reporting
# ============================================================

class TestReportInvalidNofos:
    """Tests detect_invalid_nofos → report_invalid_nofos chain."""

    def test_writes_report_csv(self, tmp_path):
        df = pd.DataFrame({
            'program_name': ['Program1', 'Program2'],
            'nofo': ['bad_nofo1', 'bad_nofo2; bad_nofo3'],
            'award': ['1R01CA123456', '2P01CA789012'],
        })
        report_path = tmp_path / 'invalid_nofos_report.csv'

        invalid_df = detect_invalid_nofos(df)
        report_invalid_nofos(invalid_df, str(report_path), printout=True)

        assert len(invalid_df) > 0
        assert 'invalid_nofo' in invalid_df.columns
        assert report_path.exists()
        assert not pd.read_csv(report_path).empty

    def test_no_file_when_all_valid(self, tmp_path, sample_dataframe):
        report_path = tmp_path / 'invalid_nofos_report.csv'

        invalid_df = detect_invalid_nofos(sample_dataframe)
        report_invalid_nofos(invalid_df, str(report_path), printout=True)

        assert invalid_df.empty
        assert not report_path.exists()


class TestReportInvalidAwards:
    """Tests detect_invalid_awards → report_invalid_awards chain."""

    def test_writes_report_csv(self, tmp_path):
        df = pd.DataFrame({
            'program_name': ['Program1', 'Program2'],
            'nofo': ['PAR-22-123', 'RFA-CA-22-456'],
            'award': ['bad_award1', 'bad_award2; bad_award3'],
        })
        report_path = tmp_path / 'invalid_awards_report.csv'

        invalid_df = detect_invalid_awards(df)
        report_invalid_awards(invalid_df, str(report_path), printout=True)

        assert len(invalid_df) > 0
        assert 'invalid_award' in invalid_df.columns
        assert report_path.exists()
        assert not pd.read_csv(report_path).empty

    def test_no_file_when_all_valid(self, tmp_path, sample_dataframe):
        report_path = tmp_path / 'invalid_awards_report.csv'

        invalid_df = detect_invalid_awards(sample_dataframe)
        report_invalid_awards(invalid_df, str(report_path), printout=True)

        assert invalid_df.empty
        assert not report_path.exists()


# ============================================================
# Column Validation and Renaming
# ============================================================

class TestValidateAndRenameColumns:

    def test_successful_rename(self):
        df = pd.DataFrame({
            'Old Column 1': [1, 2, 3],
            'Old Column 2': ['A', 'B', 'C'],
            'Old Column 3\r\n(with newline)': [4, 5, 6],
            'login_id': ['user1', 'user2', 'user3'],
        })
        col_dict = {
            'Old Column 1': 'new_col1',
            'Old Column 2': 'new_col2',
            'Old Column 3\n(with newline)': 'new_col3',
            'login_id': 'login_id',
        }

        result = validate_and_rename_columns(df, col_dict)
        assert list(result.columns) == ['new_col1', 'new_col2', 'new_col3']

    def test_column_mismatch_raises(self):
        df = pd.DataFrame({
            'Wrong Column 1': [1, 2, 3],
            'Wrong Column 2': ['A', 'B', 'C'],
        })
        col_dict = {
            'Correct Column 1': 'new_col1',
            'Correct Column 2': 'new_col2',
        }

        with pytest.raises(ValueError, match='Column names do not match expected'):
            validate_and_rename_columns(df, col_dict)

    def test_missing_login_id_last_column_raises(self):
        """Columns match col_dict, but after renaming the last column
        is not 'login_id' — so the function raises."""
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


# ============================================================
# DataFrame Cleaning Helpers
# ============================================================

class TestReplaceBlankValues:

    def test_replaces_none_and_empty(self):
        df = pd.DataFrame({
            'column1': [None, 'value', ''],
            'column2': [None, 'data', ''],
        })
        blank_dict = {'column1': 'Default1', 'column2': 'Default2'}

        result = replace_blank_values(df, blank_dict)

        assert result['column1'].tolist() == ['Default1', 'value', 'Default1']
        assert result['column2'].tolist() == ['Default2', 'data', 'Default2']


class TestForceReplaceCommaSeparation:

    def test_replaces_commas_in_listed_columns_only(self):
        df = pd.DataFrame({
            'program_name': ['Program1', 'Program, comma'],
            'focus_area': ['Area1,Area2', 'Area3, Area4'],
            'cancer_type': ['Type1,Type2', 'Type3, Type4'],
        })

        result = force_replace_comma_separation(df, ['focus_area',
                                                      'cancer_type',
                                                      'invalid_col'])

        # program_name should be untouched
        assert result['program_name'].tolist() == ['Program1', 'Program, comma']
        assert result['focus_area'].tolist() == ['Area1;Area2', 'Area3; Area4']
        assert result['cancer_type'].tolist() == ['Type1;Type2', 'Type3; Type4']


class TestRemoveExtraDefaultListValues:

    def test_removes_default_when_other_values_present(self):
        df = pd.DataFrame({
            'focus_area': ['Multiple Areas;Specific Area',
                           'Specific Area',
                           'Multiple Areas'],
            'program_acronym': ['PROG1', 'PROG2', 'PROG3'],
        })
        blank_dict = {'focus_area': 'Multiple Areas'}

        result = remove_extra_default_list_values(df, blank_dict)

        assert result['focus_area'].tolist() == [
            'Specific Area',
            'Specific Area',
            'Multiple Areas',
        ]


# ============================================================
# Program ID and Duplicate Detection
# ============================================================

class TestCreateProgramId:

    def test_various_inputs(self):
        assert create_program_id('Test Program') == 'test_program'
        assert create_program_id('Test-Program') == 'test_program'
        assert create_program_id('Test Program!') == 'test_program'
        assert create_program_id('UPPER Case') == 'upper_case'
        assert create_program_id('Non-Alpha & Num*^()Char') == 'non_alpha_num_char'


class TestGenerateProgramIdColumn:

    def test_prefers_acronym_falls_back_to_name(self):
        df = pd.DataFrame({
            'program_acronym': ['PROG1', '', 'NAN'],
            'program_name': ['First Program', 'Second Program', 'Third Program'],
        })

        result = generate_program_id_column(df)
        assert result['program_id'].tolist() == ['prog1', 'second_program',
                                                  'third_program']


class TestCheckForDuplicateNames:

    def test_duplicates_user_continues(self, monkeypatch):
        df = pd.DataFrame({
            'program_name': ['Dup Program', 'Dup Program', 'Unique'],
            'program_acronym': ['DUP', 'DUP', 'UNIQ'],
        })
        monkeypatch.setattr('builtins.input', lambda _: 'y')

        assert check_for_duplicate_names(df) is True

    def test_duplicates_user_aborts(self, monkeypatch):
        df = pd.DataFrame({
            'program_name': ['Dup Program', 'Dup Program', 'Unique'],
            'program_acronym': ['DUP', 'DUP', 'UNIQ'],
        })
        monkeypatch.setattr('builtins.input', lambda _: 'n')

        assert check_for_duplicate_names(df) is False

    def test_no_duplicates_returns_true(self):
        df = pd.DataFrame({
            'program_name': ['Unique 1', 'Unique 2'],
            'program_acronym': ['U1', 'U2'],
        })

        assert check_for_duplicate_names(df) is True


# ============================================================
# Value Fixing and User Prompts
# ============================================================

class TestApplyValueFix:

    def test_replaces_matched_values(self):
        programs_df = pd.DataFrame({
            'program_name': ['Prog A', 'Prog B'],
            'nofo': ['BAD-123', 'GOOD-456'],
        })
        reviewed_df = pd.DataFrame({
            'program_name': ['Prog A'],
            'invalid_nofo': ['BAD-123'],
            'suggested_fix': ['GOOD-123'],
        })

        result = apply_value_fix(programs_df, reviewed_df, 'nofo')

        assert result.loc[0, 'nofo'] == 'GOOD-123'
        assert result.loc[1, 'nofo'] == 'GOOD-456'


class TestPromptToContinue:

    def test_empty_df_returns_true(self):
        assert prompt_to_continue(pd.DataFrame()) is True

    def test_non_empty_user_continues(self, monkeypatch):
        df = pd.DataFrame({
            'program_name': ['Program1', 'Program2'],
            'invalid_nofo': ['bad1', 'bad2'],
        })
        monkeypatch.setattr('builtins.input', lambda _: 'Y')

        assert prompt_to_continue(df) is True

    def test_non_empty_user_aborts(self, monkeypatch):
        df = pd.DataFrame({
            'program_name': ['Program1', 'Program2'],
            'invalid_nofo': ['bad1', 'bad2'],
        })
        monkeypatch.setattr('builtins.input', lambda _: 'N')

        assert prompt_to_continue(df) is False


# ============================================================
# Integration: load_and_clean_programs
# ============================================================

class TestLoadAndCleanPrograms:

    def test_happy_path(self, valid_program_csv, sample_column_dict):
        continue_bool, result_df = load_and_clean_programs(
            valid_program_csv, sample_column_dict)

        assert continue_bool is True
        assert isinstance(result_df, pd.DataFrame)
        assert not result_df.empty

    def test_user_abort(self, valid_program_csv, sample_column_dict):
        with (
            patch('modules.gather_program_data.check_for_duplicate_names',
                  return_value=False),
            patch('modules.gather_program_data.prompt_to_continue',
                  return_value=False),
        ):
            continue_bool, result_df = load_and_clean_programs(
                valid_program_csv, sample_column_dict)

            assert continue_bool is False
            assert isinstance(result_df, pd.DataFrame)

    def test_with_nofo_and_award_corrections(self, sample_column_dict, tmp_path):
        """Loads a CSV with invalid NOFOs/Awards, applies corrections from
        reviewed CSVs, and verifies the pipeline still succeeds."""
        test_dir = os.path.join(os.path.dirname(__file__), "test_files")
        csv_path = os.path.join(test_dir, "programs_invalid_nofo_award.csv")
        nofo_corrections = os.path.join(test_dir, "programs_nofo_corrections.csv")
        award_corrections = os.path.join(test_dir, "programs_award_corrections.csv")

        with (
            patch('config.REVIEWED_NOFO_INPUT', nofo_corrections),
            patch('config.REVIEWED_AWARD_INPUT', award_corrections),
            patch('config.INVALID_NOFOS_REPORT', str(tmp_path / 'nofos.csv')),
            patch('config.INVALID_AWARD_REPORT', str(tmp_path / 'awards.csv')),
        ):
            continue_bool, result_df = load_and_clean_programs(
                csv_path, sample_column_dict)

            assert continue_bool is True
            assert isinstance(result_df, pd.DataFrame)


# ============================================================
# Integration: gather_program_data
# ============================================================

class TestGatherProgramData:

    def test_successful_output(self, tmp_path, valid_program_csv,
                               sample_column_dict):
        tmp_output = tmp_path / 'programs_intermediate.csv'

        with (
            patch('config.QUALTRICS_COLS', sample_column_dict),
            patch('config.PROGRAMS_INTERMED_PATH', tmp_output),
            patch('os.makedirs', wraps=os.makedirs),
            patch('sys.exit') as mock_exit,
        ):
            result_df = gather_program_data(valid_program_csv)
            saved_df = pd.read_csv(tmp_output)

            assert not result_df.empty
            assert tmp_output.exists()
            pd.testing.assert_frame_equal(saved_df, result_df)
            mock_exit.assert_not_called()

    def test_user_exit(self, valid_program_csv, sample_column_dict):
        with (
            patch('config.QUALTRICS_COLS', sample_column_dict),
            patch('sys.exit') as mock_exit,
            patch('modules.gather_program_data.load_and_clean_programs',
                  return_value=(False, pd.DataFrame())),
        ):
            gather_program_data(valid_program_csv)
            mock_exit.assert_called_once()

    def test_regression_against_expected_output(self, test_qualtrics_csv, tmp_path):
        """Golden-file regression test: run gather_program_data with real
        Qualtrics sample and compare output to a known-good CSV."""
        expected_csv = os.path.join(os.path.dirname(__file__),
                                    "test_files", "expected",
                                    "test_program_gathering_output.csv")
        full_column_dict = {
            "Name of Key Program": "program_name",
            "Acronym for key program": "program_acronym",
            "Focus Area (select all that apply)": "focus_area",
            "DOC": "doc",
            "Primary Contact (PI)": "contact_pi",
            "Primary Contact (PI) email": "contact_pi_email",
            "NIH Contact (Program Officer/Program Director)": "contact_nih",
            "NIH Contact (Program Officer/Program Director) email":
                "contact_nih_email",
            'NOFO number (eg. format as "RFA-CA-00-000") '
            '(If more than one, separate with ; semicolon)': "nofo",
            "Grant/Award number {parent award FORMAT LL#CA######, "
            "eg. UG3CA260607} (If more than one, separate with ; semicolon)":
                "award",
            "Link to program website": "program_link",
            "Link to data or DCC if available": "data_link",
            "What type of cancer is the primary focus of the program? "
            "(Check all that\napply)": "cancer_type",
            "Login ID": "login_id",
        }
        tmp_output = tmp_path / 'test_actual_program_output.csv'

        with (
            patch('config.QUALTRICS_COLS', full_column_dict),
            patch('config.PROGRAMS_INTERMED_PATH', tmp_output),
            patch('config.INVALID_AWARD_REPORT',
                  str(tmp_path / 'invalid_awards.csv')),
            patch('modules.gather_program_data.prompt_to_continue',
                  return_value=True),
        ):
            result_df = gather_program_data(test_qualtrics_csv)
            expected_df = pd.read_csv(expected_csv)

            assert not result_df.empty
            pd.testing.assert_frame_equal(result_df, expected_df)
            assert tmp_output.exists()
