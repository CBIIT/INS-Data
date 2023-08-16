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

# TODO: Add detection and resolution for bad NOFO lists (e.g. ;; or PA18-:91)

import pandas as pd


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
    assert False, f"Key value '{key_value}' not found in file."


def load_and_clean_programs(csv_filepath: str, col_dict: dict) -> pd.DataFrame:
    """Load and clean Key Programs data from a Qualtrics CSV.

    :param csv_filepath: Path to the CSV file to load. Defined in config.py
    :type csv_filepath: str

    :param col_dict: Dictionary of column names to expect within CSV and new
                    names to assign to each. {old_name: new_name, ...}
    :type col_dict: dict

    :return: The cleaned Key Programs data as a pandas DataFrame.
    :rtype: pd.DataFrame
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
    # TODO: Consider moving assertion to formal tests
    actual_cols = df.columns.tolist()
    expected_cols = list(col_dict.keys())
    assert (actual_cols == expected_cols), (
        f"Column names do not match expected."
        f"Expected columns: \n {expected_cols} \n ---"
        f"Actual columns: \n {actual_cols}")


    # Rename columns with defined dictionary
    df = df.rename(columns=col_dict)

    # Drop second header row with survey question IDs
    df = df.drop(axis=0, index=0).reset_index(drop=True)

    # Check that the last column is login_id and then drop it
    # TODO: Again, consider moving assertion to formal tests
    assert (df.columns[-1] == "login_id"), ( 
            f"Unexpected final column: {df.columns[-1]}")
    df = df.drop(columns=["login_id"])

    return df


