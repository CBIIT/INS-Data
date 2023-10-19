# data_preparation.py
# 2023-07-27 ZD 
#
# This script defines primary function `load_and_clean_programs` that will 
# input a versioned `qualtrics_output_{YYYY-MM-DD}.csv` file, perform cleaning 
# and processing steps, and then output a `key_programs` pandas DataFrame (df).
#
# It also defines helper function `find_header_location`, which is called 
# to identify where the desired data begins within the CSV. i.e. it identifies
# the location of the desired "top-left, A1 cell" and extracts the data of 
# interest from that anchor location.
#
# The output `key_programs` df contains the following columns: 
# 'program_name'
# 'program_acronym'
# 'focus_area'
# 'doc'
# 'contact_pi'
# 'contact_pi_email'
# 'contact_nih'
# 'contact_nih_email'
# 'nofo'
# 'award'
# 'program_link'
# 'data_link'
# 'cancer_type

import pandas as pd
import re
import os
import sys
# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config



def find_header_location(csv_filepath:str, key_value:str) -> 'tuple[int,int]':
    """Detect the row and column where the given key_value is found.

    :param csv_filepath: Path to the CSV file to load. Defined in config.py.
    :type csv_filepath: str

    :param key_value: String value to search for within CSV.
    :type key_value: str

    :return: The row, column location of the key value within the CSV.
    :rtype: tuple[int,int]
    """

    with open(csv_filepath, 'r') as file:
        for row, line in enumerate(file):
            # Read the first row (header row) and split it by commas
            header_row = line.strip().split(',')  
            for col, header_value in enumerate(header_row):
                if key_value == header_value.strip():
                    return row, col

    # If the loop finishes without finding the key_value, raise an Error
    raise ValueError(f"Key value '{key_value}' not found in file.")




def clean_nofo_and_award_cols(df:pd.DataFrame) -> pd.DataFrame:
    """Clean obvious formatting mistakes from NOFO and Award strings."""

    # Ensure 'nofo' and 'award' columns are present in the DataFrame
    if 'nofo' not in df.columns or 'award' not in df.columns:
        raise ValueError(f"Columns 'nofo' and 'award' must be present in "
                         f"the DataFrame.")

    # Clean 'nofo' column. Remove spaces and replace any commas with semicolon
    df['nofo'] = df['nofo'].apply(lambda x: str(x).replace(" ", "")
                                  .replace(",", ";") if pd.notna(x) else x)

    # Clean 'award' column. Remove spaces and replace any commas with semicolon
    df['award'] = df['award'].apply(lambda x: str(x).replace(" ", "")
                                    .replace(",", ";") if pd.notna(x) else x)

    return df



def clean_and_split_nofo_award_strings(value_list):
    """Parse string list into actual [list], separating by semicolon"""

    if pd.isna(value_list) or value_list.lower() == 'nan':
        return []
    return [nf.strip() for nf in str(value_list).split(';') if nf.strip()]



def is_valid_nofo(nofo:str) -> bool:
    """Check if the provided NOFO follows a valid format.

    Valid NOFO Formats:
    - RFA
        - RFA-XX-##-###
    - PAR / OTA
        - PAR-##-###
        - PAR##-###
        - OTA##-###
        - OTA-##-###
    - No Prefix
        - XX##-###
        - XX-##-###

    Where "X" is any letter, "#" is any number, "-" is a hyphen, and any explicit
    letters ("PAR", "OTA", "RFA") must be those exact letters.

    :param str nofo: NOFO string to be validated.
    :return bool: True if the NOFO is valid, False otherwise.
    """

    # Define the valid NOFO regex pattern
    pattern_no_prefix = re.compile(r'^[A-Z]{2}-?\d{2}-\d{3}$')
    pattern_rfa = re.compile(r'^RFA-[A-Z]{2}-\d{2}-\d{3}$')
    pattern_par_ota = re.compile(r'^(PAR|OTA)-?\d{2}-\d{3}$')

    # Check if the provided nofo matches the valid pattern
    return any((pattern_no_prefix.match(nofo),
                pattern_rfa.match(nofo),
                pattern_par_ota.match(nofo)))



def validate_nofos(program_name:str, nofo_list:list) -> pd.DataFrame:
    """Validate that provided NOFOs fit the expected regex format. Output any 
    potential errors for manual review and resolution.

    :param str program_name: Name of associated program
    :param list nofo_list: Parsed list of NOFO strings
    :return pd.DataFrame: Pandas DataFrame of invalid NOFOs and programs, or empty
    """
    # Create list of invalid NOFO strings
    invalid_nofos = [nf for nf in nofo_list if not is_valid_nofo(nf) 
                                                and nf.lower() != 'nan']

    # Create dataframe of any invalid NOFOs and associated program
    if invalid_nofos:
        invalid_df = pd.DataFrame({'program_name': program_name, 
                                   'invalid_nofo': invalid_nofos})
        return invalid_df
    
    # Return empty dataframe if no invalid NOFOs detected
    return pd.DataFrame()



def detect_invalid_nofos(df:pd.DataFrame) -> pd.DataFrame:
    """Detect rows with potentially invalid NOFOs for manual review.

    :param pd.DataFrame df: Input DataFrame containing 'program_name' 
                            and 'nofo' columns.
    """
    # Create a new column of parsed nofos as list
    df['parsed_nofos'] = df['nofo'].apply(clean_and_split_nofo_award_strings)

    # Concatenate DataFrames row-wise
    invalid_nofos_df = pd.concat(df.apply(lambda row: validate_nofos(
                                    row['program_name'], row['parsed_nofos']),
                                    axis=1).tolist())

    return invalid_nofos_df



def report_review_continue(invalid_nofos_df:pd.DataFrame) -> bool:
    """Prompt the user to review potential issues and continue or stop the 
    workflow. If no issues are found, then skip the prompt. 

    :param pd.DataFrame invalid_nofos_df: DataFrame containing invalid NOFOs.
    :return bool: True if the user wants to continue, False otherwise.
    """
    if not invalid_nofos_df.empty:
        print("Potentially invalid NOFOs found:")
        print(invalid_nofos_df)

        # Export invalid NOFOs to csv
        validity_reports_dir = config.REPORTS_DIR
        if not os.path.exists(config.REPORTS_DIR):
            os.makedirs(validity_reports_dir)
        invalid_nofos_path = os.path.join(validity_reports_dir, 
                                          config.INVALID_NOFOS_FILENAME)
        invalid_nofos_df.to_csv(invalid_nofos_path, index=False)

        # Prompt the user to review potential issues
        review_input = input(f"---\n\n"
                            f"Potentially invalid NOFO report saved to: "
                            f"{invalid_nofos_path}\n"
                            f"Please review potential issues and consider "
                            f"manual fixes to the Qualtrics CSV. \n"
                            f"CONTINUE ANYWAY? (Y/N): ").upper()
        # If user inputs 'Y', then return True
        return review_input == 'Y'
    else:
        print("All good! No potentially invalid NOFOs detected.")
        return True



def load_and_clean_programs(csv_filepath: str, col_dict: dict) -> (bool, pd.DataFrame):
    """Load and clean Key Programs data from a Qualtrics CSV.

    :param str csv_filepath: Path to the CSV file to load. Defined in config.py
    :param dict col_dict: Dictionary of column names to expect within CSV and new
                    names to assign to each. {old_name: new_name, ...}
    :return bool continue_bool: True/False indicator to continue or stop process
    :return pd.DataFrame df: The cleaned Key Programs data as a pandas DataFrame.
    """

    # Get string of first key column to look for within file
    # e.g. "Name of Key Program"
    key_value = list(col_dict.keys())[0]

    # Get row,col of given header text
    # This will be used to determine how many irrelevant rows and columns to 
    # skip/drop when processing the CSV.
    header_row, header_col = find_header_location(csv_filepath, key_value)

    # Load file, skip leading rows, and then drop leading columns
    df = (pd.read_csv(csv_filepath, skiprows=header_row)
      .iloc[:, header_col:])
    
    # Use dictionary to check for unexpected or missing columns in data
    actual_cols = df.columns.tolist()
    expected_cols = list(col_dict.keys())
    if actual_cols != expected_cols:
        raise ValueError(
            f"Column names do not match expected.\n"
            f"Expected columns: \n {expected_cols} \n"
            f"Actual columns: \n {actual_cols}")

    # Rename columns with defined dictionary
    df = df.rename(columns=col_dict)

    # Drop second header row with survey question IDs
    df = df.drop(axis=0, index=0).reset_index(drop=True)

    # Remove spaces and replace commas in NOFO and Award columns
    df = clean_nofo_and_award_cols(df)

    # Remove leading/trailing whitespace from program names
    df['program_name'] = df['program_name'].str.strip()

    # Detect potentially invalid NOFO values
    print("Checking for NOFO validity...")

    # Make a copy to prevent changing main df
    nofo_validity_df = df.copy()
    invalid_nofos_df = detect_invalid_nofos(nofo_validity_df)

    # Prompt user to review and stop or continue if issues are found
    continue_bool = report_review_continue(invalid_nofos_df)
        
    # Check that the last column is login_id and then drop it
    if df.columns[-1] != "login_id":
        raise ValueError(f"Unexpected final column: {df.columns[-1]}")
    df = df.drop(columns=["login_id"])

    return continue_bool, df



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # STEP 1: Data Prep - Load and clean Key Programs data
    print(f"Loading and processing {config.QUALTRICS_CSV_PATH}...")

    # Load and clean Key Programs
    continue_bool, key_programs_df = load_and_clean_programs(
                                    csv_filepath = config.QUALTRICS_CSV_PATH, 
                                    col_dict = config.QUALTRICS_COLS)

    if continue_bool == True:
        # Export cleaned Key Programs file
        key_programs_df.to_csv(config.CLEANED_KEY_PROGRAMS_CSV, index=False)
        print(f"Success! Saved {config.CLEANED_KEY_PROGRAMS_CSV}.")

    else:
        sys.exit("Process ended by user. Cleaned Key Programs csv not saved. "
            "Make manual edits to input CSV and retry.")
