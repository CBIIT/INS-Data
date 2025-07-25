"""
gather_program_data.py
2023-07-27 ZD 

This script defines primary function `gather_program_data` that will 
input a versioned qualtrics CSV containing curated program data. It will read 
configured values, perform validation, cleaning, and processing steps, and 
then output a DataFrame and matching CSV file containing programs.

The output dataframe contains the following columns: 
    program_id
    program_name
    program_acronym
    focus_area
    doc
    contact_pi
    contact_pi_email
    contact_nih
    contact_nih_email
    nofo,award
    program_link
    data_link
    cancer_type
"""

import os
import re
import sys
import pandas as pd
# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config



def find_header_location(csv_filepath:str, key_value:str) -> 'tuple[int,int]':
    """Detect the row and column where the given key_value is found.

        Args:
            csv_filepath (str): Path to the CSV file to load. Defined in 
                config.py.
            key_value (str): String value to search for within CSV.

        Returns:
            tuple[int, int]: The row, column location of the key value within 
                the CSV.
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



def drop_obsolete_columns(df:pd.DataFrame, obsolete_str:str="obsolete"):
    """Drop columns with header containing 'obsolete' string added during 
    config renaming."""

    # Identify any column headers renamed to include 'obsolete' string
    obsolete_cols = [col for col in df.columns if obsolete_str in col]

    # If any are found, drop them from the dataframe
    if obsolete_cols:
        df = df.drop(obsolete_cols, axis=1)
    
    return df



def clean_nofo_and_award_cols(df:pd.DataFrame) -> pd.DataFrame:
    """Clean obvious formatting mistakes from NOFO and Award strings."""

    # Ensure 'nofo' and 'award' columns are present in the DataFrame
    if 'nofo' not in df.columns or 'award' not in df.columns:
        raise ValueError(f"Columns 'nofo' and 'award' must be present in "
                         f"the DataFrame.")

    # Clean 'nofo' column. Remove spaces and replace any commas with semicolon
    df['nofo'] = df['nofo'].apply(lambda x: str(x)
                                  .replace(" ", "")
                                  .replace("\u00A0", "") # Non-breaking space
                                  .replace(",", ";")

                                  if pd.notna(x) else x)

    # Clean 'award' column. Remove spaces and replace any commas with semicolon
    df['award'] = df['award'].apply(lambda x: str(x)
                                    .replace(" ", "")
                                    .replace("\u00A0", "") # Non-breaking space
                                    .replace(",", ";")
                                    if pd.notna(x) else x)

    return df



def clean_and_split_nofo_award_strings(value_list):
    """Parse string list into an actual [list], separating by a semicolon."""

    # Check if value_list is empty or contains 'nan'
    if pd.isna(value_list) or value_list.lower() == 'nan':
        return []

    # Split string by semicolon, strip whitespaces, filter out empty strings
    cleaned_values = [nf.strip() for nf in str(value_list)
                      .split(';') if nf.strip()]

    return cleaned_values



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

    Args:
        nofo (str): NOFO string to be validated.

    Returns:
        bool: True if the NOFO is valid, False otherwise.
    """

    # Define the valid NOFO regex pattern
    pattern_no_prefix = re.compile(r'^[A-Z]{2}-?\d{2}-\d{3}$')
    pattern_rfa = re.compile(r'^RFA-[A-Z]{2}-\d{2}-\d{3}$')
    pattern_par_ota = re.compile(r'^(PAR|OTA)-?\d{2}-\d{3}$')

    # Check if the provided nofo matches the valid pattern
    return any((pattern_no_prefix.match(nofo),
                pattern_rfa.match(nofo),
                pattern_par_ota.match(nofo)))



def is_valid_award(award: str) -> bool:
    """Check if the provided award follows a valid format.

    Valid Award Formats:
        - Grants
            1.   X##XX######-##A#S#   Supplement (e.g. 01A1S1)
            2.   X##XX######-##X#     Supplement (e.g. 01S1 or 01A1)
            3.   X##XX######-##       Yearly (e.g. 01)
            4.   X##XX######
            5.   XX######
        - Contracts
            - Contract formats vary greatly, so a loose pattern is used
            - Must start with 6 consecutive numbers
            - Must be at least 16 characters long
            - Must contain only numbers, letters, and hyphens
          
    Where X is any letter, # is any number, and any other character must be an 
    explicit match ("-", "A", "S"). All grant types can be optionally preceeded
    by a single number digit and still be valid. 

    Args:
        award (str): Award string to be validated.

    Returns:
        bool: True if the award is valid, False otherwise.
    """

    # Define the valid award regex patterns
    pattern_grant_type1 = re.compile(r'^[A-Z]\d{2}[A-Z]{2}\d{6}-\d{2}[A]\d[S]\d$')
    pattern_grant_type1a = re.compile(r'^\d{1}[A-Z]\d{2}[A-Z]{2}\d{6}-\d{2}[A]\d[S]\d$')
    pattern_grant_type2 = re.compile(r'^[A-Z]\d{2}[A-Z]{2}\d{6}-\d{2}[A-Z]\d$')
    pattern_grant_type2a = re.compile(r'^\d{1}[A-Z]\d{2}[A-Z]{2}\d{6}-\d{2}[A-Z]\d$')
    pattern_grant_type3 = re.compile(r'^[A-Z]\d{2}[A-Z]{2}\d{6}-\d{2}$')
    pattern_grant_type3a = re.compile(r'^\d{1}[A-Z]\d{2}[A-Z]{2}\d{6}-\d{2}$')
    pattern_grant_type4 = re.compile(r'^[A-Z]\d{2}[A-Z]{2}\d{6}$')
    pattern_grant_type4a = re.compile(r'^\d{1}[A-Z]\d{2}[A-Z]{2}\d{6}$')
    pattern_grant_type5 = re.compile(r'^[A-Z]{2}\d{6}$')
    pattern_grant_type5a = re.compile(r'^\d{1}[A-Z]{2}\d{6}$')
    pattern_contract_type = re.compile(r'^\d{6,}[A-Za-z0-9-]{10,}$')

    # Check if the provided award matches any of the valid patterns
    return any((pattern_grant_type1.match(award),
                pattern_grant_type1a.match(award),
                pattern_grant_type2.match(award),
                pattern_grant_type2a.match(award),
                pattern_grant_type3.match(award),
                pattern_grant_type3a.match(award),
                pattern_grant_type4.match(award),
                pattern_grant_type4a.match(award),
                pattern_grant_type5.match(award),
                pattern_grant_type5a.match(award),
                pattern_contract_type.match(award)))



def validate_nofos(program_name:str, nofo_list:list) -> pd.DataFrame:
    """Validate that provided NOFOs fit the expected regex format. Output any
    potential issues for manual review and correction.

    Args:
        program_name (str): Name of the associated program.
        nofo_list (list): Parsed list of NOFO strings.

    Returns:
        pd.DataFrame: DataFrame of invalid NOFOs and programs, or empty.
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



def validate_awards(program_name: str, award_list: list) -> pd.DataFrame:
    """Validate that provided awards fit the expected regex format. Output any
    potential issues for manual review and correction.

    Args:
        program_name (str): Name of the associated program.
        award_list (list): Parsed list of award strings.

    Returns:
        pd.DataFrame: DataFrame of invalid awards and programs, or empty.
    """
    # Create a list of invalid award strings
    invalid_awards = [award for award in award_list if not is_valid_award(award)
                      and award.lower() != 'nan']

    # Create a DataFrame of any invalid awards and associated programs
    if invalid_awards:
        invalid_df = pd.DataFrame({'program_name': program_name,
                                   'invalid_award': invalid_awards})
        return invalid_df
    
    # Return an empty DataFrame if no invalid awards are detected
    return pd.DataFrame()



def detect_invalid_nofos(df:pd.DataFrame) -> pd.DataFrame:
    """Detect rows with potentially invalid NOFOs for manual review.

    Args:
        df (pd.DataFrame): Input DataFrame containing 'program_name' and 
            'nofo' columns.

    Returns:
        pd.DataFrame: DataFrame of potentially invalid NOFOs and programs.
    """

    # Create a new column of parsed nofos as list
    df['parsed_nofos'] = df['nofo'].apply(clean_and_split_nofo_award_strings)

    # Concatenate DataFrames row-wise
    invalid_nofos_df = pd.concat(df.apply(lambda row: validate_nofos(
                                    row['program_name'], row['parsed_nofos']),
                                    axis=1).tolist())

    return invalid_nofos_df
    


def detect_invalid_awards(df: pd.DataFrame) -> pd.DataFrame:
    """Detect rows with potentially invalid awards for manual review.

    Args:
        df (pd.DataFrame): Input DataFrame containing 'program_name' and 
            'award' columns.

    Returns:
        pd.DataFrame: DataFrame of potentially invalid awards and programs.
    """

    # Create a new column of parsed awards as a list
    df['parsed_awards'] = df['award'].apply(clean_and_split_nofo_award_strings)

    # Concatenate DataFrames row-wise
    invalid_awards_df = pd.concat(df.apply(lambda row: validate_awards(
                                    row['program_name'], row['parsed_awards']),
                                    axis=1).tolist())

    return invalid_awards_df



def report_invalid_nofos(invalid_nofos_df:pd.DataFrame, 
                         report_path:str,
                         printout:bool = False) -> None:
    """Print and export output of invalid NOFOs found within key programs df.

    Args:
        invalid_nofos_df (pd.DataFrame): DataFrame containing invalid NOFOs.
        report_path (str): Filepath for CSV export of the invalid NOFO report.
        printout (bool): If true, print output to the terminal. Default False.
    """

    if not invalid_nofos_df.empty:
        print(f"{len(invalid_nofos_df)} potentially invalid NOFOs found.")
        if printout:
            print(invalid_nofos_df)

        # Export invalid NOFOs to csv
        os.makedirs(os.path.dirname(report_path), exist_ok=True)
        invalid_nofos_df.to_csv(report_path, index=False)
        print(f"Invalid NOFOs saved for review "
              f"and correction in {report_path}.")
    else: 
        print("All good. No potentialy invalid NOFOs detected.")
    
    return None



def report_invalid_awards(invalid_awards_df: pd.DataFrame,
                          report_path: str,
                          printout: bool = False) -> None:
    """Print and export the output of invalid awards found within key programs df.

    Args:
        invalid_awards_df (pd.DataFrame): DataFrame containing invalid awards.
        report_path (str): Filepath for CSV export of the invalid awards report.
        printout (bool): If true, print output to the terminal. Default False.
    """

    if not invalid_awards_df.empty:
        print(f"{len(invalid_awards_df)} potentially invalid awards found.")
        if printout:
            print(invalid_awards_df)

        # Export invalid awards to csv
        os.makedirs(os.path.dirname(report_path), exist_ok=True)
        invalid_awards_df.to_csv(report_path, index=False)
        print(f"Invalid awards saved for review "
              f"and correction in {report_path}.")
    else: 
        print("All good. No potentially invalid awards detected.")
    
    return None



def prompt_to_continue(invalid_df: pd.DataFrame) -> bool:
    """Provide the user with a prompt to continue or stop the workflow if 
        issues are present.

    Args:
        invalid_df (pd.DataFrame): DataFrame containing invalid data.

    Returns:
        bool: True if the user wants to continue, False otherwise.
    """

    # If no invalid values detected, skip the user prompt
    if invalid_df.empty:
        return True
    
    else:
        # Parse award or nofo type from "invalid_{col_name}"
        val_type = invalid_df.columns[1].split('_')[-1].capitalize()

        print(f"---\n"
            f"Please review {len(invalid_df)} potential "
            f"{val_type} issues. \n"
            f"Consider manual fixes to the Qualtrics CSV or creating a "
            f"versioned invalid{val_type}Report_reviewed.csv` "
            f"in the data/reviewed/ directory")
        
        continue_bool = input(f"\n\tContinue with known {val_type} "
                              f"issues? (Y/N): "
                               ).upper()
        return continue_bool == 'Y'
    


def apply_value_fix(key_programs_df: pd.DataFrame, 
                    reviewed_df: pd.DataFrame, 
                    column_name: str) -> pd.DataFrame:
    """Apply fixes from the reviewed file to the key programs DataFrame.

    Args:
        key_programs_df (pd.DataFrame): The original key programs DataFrame.
        reviewed_df (pd.DataFrame): The DataFrame containing fixes.
        column_name (str): The name of the column to apply fixes to.

    Returns:
        pd.DataFrame: The key programs DataFrame with applied fixes.
    """

    # Iterate through each row in the reviewed DataFrame
    for index, row in reviewed_df.iterrows():
        program_name = row['program_name']
        invalid_value = row[f'invalid_{column_name}']
        suggested_fix = row['suggested_fix']

        # Find the row in key_programs_df where program_name matches
        program_rows = key_programs_df[
                        key_programs_df['program_name'] == program_name]

        # Iterate through rows to update invalid_value with suggested_fix
        for idx, program_row in program_rows.iterrows():
            # Check if the invalid_value is present in the current row
            if invalid_value in program_row[column_name]:
                # Apply the suggested_fix
                key_programs_df.at[idx, column_name] = program_row[
                    column_name].replace(invalid_value, suggested_fix)

    return key_programs_df



def validate_and_rename_columns(df: pd.DataFrame, 
                                col_dict: dict) -> pd.DataFrame:
    """Validate and rename columns based on col_dict.

    Args:
        df (pd.DataFrame): Input DataFrame.
        col_dict (dict): Dictionary of column names to expect within the CSV 
            and new names to assign to each. {old_name: new_name, ...}

    Returns:
        pd.DataFrame: DataFrame with validated and renamed columns.
    """

    # Replace any carriage returns newlines with newlines in headers
    df.columns = df.columns.str.replace('\r\n', '\n')

    # Validate that columns received match expected
    actual_cols = df.columns.tolist()
    expected_cols = list(col_dict.keys())
    if actual_cols != expected_cols:
        raise ValueError(f"Column names do not match expected.\n"
                         f"Expected columns: {expected_cols}\n"
                         f"Actual columns: {actual_cols}")
    
    # Rename columns based on dictionary
    df.rename(columns=col_dict, inplace=True)

    # Check that the last column is login_id and then drop it
    if df.columns[-1] != "login_id":
        raise ValueError(f"Unexpected final column: {df.columns[-1]}")
    df = df.drop(columns=["login_id"])

    return df



def replace_blank_values(df: pd.DataFrame, blank_filler_dict: dict):
    """Replace blank values within specific columns using config file."""

    # Iterate through renaming dictionary defined in config
    for column, filler in blank_filler_dict.items():
        # Fill blank and NaN values
        df[column] = df[column].fillna(filler)
        df[column] = df[column].replace('',filler)

    return df



def force_replace_comma_separation(df: pd.DataFrame, 
                                   replace_cols: list) -> pd.DataFrame:
    """Replace any commas within specified columns with semicolons. This is
    risky and should be reviewed for unexpected splitting of list values. 

    Args:
        df (pd.DataFrame): Input DataFrame.
        replace_cols (list): List of column names to replace commas in.

    Returns:
        pd.DataFrame: DataFrame with replaced commas in specified columns.
    """

    for col in replace_cols:
        # Check if the column exists in the DataFrame
        if col in df.columns:
            # Replace commas with semicolons in the specified column
            df[col] = df[col].str.replace(',', ';')
        else:
            print(f"List-like column '{col}' not found in the DataFrame.")

    return df



def check_for_duplicate_names(df: pd.DataFrame) -> bool:
    """Check for duplicate values in program_name or program_acronym.

    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        bool: True if no duplicates are found. If duplicates found, user has
            the option to continue (True) or end (False) the process with input.
    """

    # Find rows with duplicate names or acronyms
    dup_names = df[df.duplicated(subset=['program_name'], keep=False)]
    dup_acronyms = df[df.duplicated(subset=['program_acronym'], keep=False)]

    # Check if duplicates exist before proceeding
    if not dup_names.empty or not dup_acronyms.empty:
        print(f"---\nWarning: "
              f"Duplicate values found in program_name or program_acronym.")

        # Print out any duplicates
        if not dup_names.empty:
            print("\nProgram Names with Duplicates:")
            print(dup_names[['program_name', 'program_acronym']])
        if not dup_acronyms.empty:
            print("\nProgram Acronyms with Duplicates:")
            print(dup_acronyms[['program_name', 'program_acronym']])

        print(f"\nConsider manual fixes to the Qualtrics CSV.")

        # Prompt user to continue or stop with duplicates
        continue_bool = input(f"\n\tContinue with duplicates? (Y/N): "
                              ).lower() == 'y'
        
        return continue_bool

    # If there are no duplicates, return True to continue
    return True



def create_program_id(value):
    """Create a valid program_id based on the specified column value.

    Args:
        value (str): The value from which to derive the program_id.

    Returns:
        str: The generated program_id.
    """

    # Replace hyphens and other non-alphanumeric characters with underscores
    cleaned_value = ''.join('_' if not char.isalnum() 
                            else char for char in str(value))

    # Replace consecutive non-alphanumeric characters with a single underscore
    cleaned_value = '_'.join(part for part in cleaned_value.split('_') if part)

    # Convert to lowercase
    return cleaned_value.lower()



def generate_program_id_column(df: pd.DataFrame) -> pd.DataFrame:
    """Generate a new program_id column based on the program_acronym or name.

    Args:
        df (pd.DataFrame): Input DataFrame with program_acronym and program_name
            columns.

    Returns:
        pd.DataFrame: DataFrame with the new 'program_id' column added.
    """

    # Create empty program_id column
    df['program_id'] = ''

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        program_acronym = str(row['program_acronym']).strip()

        # Check if program_acronym exists and is not a blank space
        if program_acronym and program_acronym.lower() != 'nan':
            # Use program_acronym to generate program_id
            df.at[index, 'program_id'] = create_program_id(program_acronym)
        else:
            # If no acronym, use program_name to generate program_id
            df.at[index, 'program_id'] = create_program_id(row['program_name'])

            # Print warning for missing acronym
            if pd.isna(program_acronym) or program_acronym == '':
                print(f"---\nNote: No program_acronym provided for "
                      f"Program: {row['program_name']}")

    # Move the program_id column to the first position
    cols = ['program_id'] + [col for col in df.columns if col != 'program_id']
    df = df[cols]

    return df



def remove_extra_default_list_values(df: pd.DataFrame, blank_filler_dict: dict):
    """Remove default values from lists when those values should be mutually
    exclusive with any other values. E.g. 'Multiple Cancer Types' should be 
    removed when other cancer types are also listed. 
    Args:
        df (pd.DataFrame): Input pandas dataframe
        blank_filler_dict (dict): Dictionary containing column names and default
        values. Reusing from similar `replace_blank_values` function.
    """

    # Iterate through rows in df
    for index, row in df.iterrows():

        # Iterate through dictionary of columns and values from config
        for col_name, default_value in blank_filler_dict.items():
            # Separate strings by semicolons to get values
            values = row[col_name].split(';')

            # Remove mutually exclusive default value if present with others
            if default_value in values and len(values) > 1:
                print(f"Removed extra {default_value} in {row['program_acronym']}")
                values.remove(default_value)

            # Rejoin remaining values into semicolon separated string
            df.loc[index, col_name] = ';'.join(values)

    return df



def load_and_clean_programs(csv_filepath: str, col_dict: dict) -> tuple[bool, pd.DataFrame]:
    """Load and clean Key Programs data from a Qualtrics CSV.

    Args:
        csv_filepath (str): Path to the CSV file to load. Defined in config.py
        col_dict (dict): Dictionary of column names to expect within CSV and new
                        names to assign to each. {old_name: new_name, ...}

    Returns:
        bool: True/False indicator to continue or stop the process
        pd.DataFrame: The cleaned Key Programs data as a pandas DataFrame.
    """

    # Get string of first key column to look for within file and get row
    key_value = list(col_dict.keys())[0]
    header_row, header_col = find_header_location(csv_filepath, key_value)

    # Load file, skip leading rows, and then drop leading columns
    df = pd.read_csv(csv_filepath, skiprows=header_row).iloc[:, header_col:]

    # Check and rename column names
    df = validate_and_rename_columns(df, col_dict)

    # Drop obsolete columns specified in config.py
    df = drop_obsolete_columns(df, obsolete_str="obsolete")

    # Replace specific string values with others defined in config
    for old_value, new_value in config.PROGRAM_VALUE_REPLACEMENTS.items():
        # Use regex=True to enable regular expression replacement
        df = df.replace({r'\b{}\b'.format(re.escape(old_value)): new_value}, regex=True)

    # Replace blank values in specific columns with filler values
    df = replace_blank_values(df, config.PROGRAM_BLANK_REPLACEMENTS)

    # Remove spaces and replace commas in NOFO and Award columns
    df = clean_nofo_and_award_cols(df)

    # Brute force replacement of commas with semicolons in list-like cols
    df = force_replace_comma_separation(df, ['focus_area', 'cancer_type', 'doc'])

    # Remove mutually exclusive default values if present in list-like cols
    df = remove_extra_default_list_values(df, config.PROGRAM_BLANK_REPLACEMENTS)

    # Remove leading/trailing whitespace from program names
    df['program_name'] = df['program_name'].str.strip()

    # Check for duplicate program names or acronyms
    continue_bool_names = check_for_duplicate_names(df)
    if not continue_bool_names:
        return False, df

    # Add a program_id column to the first position
    df = generate_program_id_column(df)

    # Detect potentially invalid NOFO values
    print("---\nChecking for NOFO validity...")

    # Make a copy to prevent changing main df
    nofo_validity_df = df.copy()
    invalid_nofos_df = detect_invalid_nofos(nofo_validity_df)

    # Report invalid NOFOs
    report_invalid_nofos(invalid_nofos_df, 
                         config.INVALID_NOFOS_REPORT,
                         printout=False)

    # Check if corrections to invalid NOFOs already exist
    if os.path.exists(config.REVIEWED_NOFO_INPUT):
        print(f"\nReviewed NOFO file found: {config.REVIEWED_NOFO_INPUT} \n"
              f"Applying corrections and rerunning NOFO validation...")
        
        # Read csv with NOFO corrections and then apply
        reviewed_df = pd.read_csv(config.REVIEWED_NOFO_INPUT)
        df = apply_value_fix(df, reviewed_df, 'nofo')

         # Make a copy to prevent changing main df
        nofo_validity_df = df.copy()
        invalid_nofos_df = detect_invalid_nofos(nofo_validity_df)

        # Report invalid NOFOs
        report_invalid_nofos(invalid_nofos_df, 
                             config.CORRECTED_INVALID_NOFOS_REPORT, 
                             printout=False)
        
    continue_bool_nofos = prompt_to_continue(invalid_nofos_df)

    # Detect potentially invalid 'Award' values
    print("---\nChecking for Award validity...")

    # Make a copy to prevent changing the main df
    award_validity_df = df.copy()
    invalid_awards_df = detect_invalid_awards(award_validity_df)

    # Report invalid Awards
    report_invalid_awards(invalid_awards_df, 
                          config.INVALID_AWARD_REPORT, 
                          printout=False)

    # Check if corrections to invalid Awards already exist
    if os.path.exists(config.REVIEWED_AWARD_INPUT):
        print(f"\nReviewed 'Award' file found: {config.REVIEWED_AWARD_INPUT} \n"
              f"Applying corrections and rerunning Award validation...")

        # Read CSV with 'Award' corrections and then apply
        reviewed_awards_df = pd.read_csv(config.REVIEWED_AWARD_INPUT)
        df = apply_value_fix(df, reviewed_awards_df, 'award')

        # Make a copy to prevent changing the main df
        award_validity_df = df.copy()
        invalid_awards_df = detect_invalid_awards(award_validity_df)

        # Report invalid Awards
        report_invalid_awards(invalid_awards_df, 
                              config.CORRECTED_INVALID_AWARD_REPORT, 
                              printout=False)
    continue_bool_awards = prompt_to_continue(invalid_awards_df)

    # Prompt user to continue if issues found
    continue_bool = (continue_bool_nofos and continue_bool_awards)

    return continue_bool, df



def gather_program_data(qualtrics_csv: str) -> pd.DataFrame:
    """Process a curated Qualtrics CSV containing NCI programs and associated 
    funding values. Validate, clean, and prepare for downstream data gathering.

    Args:
        qualtrics_csv (str): Filepath to the Qualtrics CSV to ingest and process

    Returns: 
        pd.DataFrame: Formatted DataFrame of programs and details
    """

    print(f"\n---\nPROGRAMS:\n"
          f"Gathering, cleaning, and saving programs data...\n---\n")

    # Load and clean programs data
    print(f"Loading and processing {qualtrics_csv}...")
    continue_bool, programs_df = load_and_clean_programs(
                                    csv_filepath = qualtrics_csv, 
                                    col_dict = config.QUALTRICS_COLS)

    if continue_bool == True:
        # Export cleaned Key Programs file
        program_filepath = config.PROGRAMS_INTERMED_PATH
        os.makedirs(os.path.dirname(program_filepath), exist_ok=True)
        programs_df.to_csv(program_filepath, index=False)
        print(f"Success! {len(programs_df)} Programs saved to {program_filepath}.")

    else:
        sys.exit(f"\n---\n"
                 f"Process ended by user. Cleaned programs CSV not saved. "
                 f"Make manual edits to input CSV and retry.")
        
    return programs_df



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Gather program data. Load, clean, and validate from curated file.
    programs_df = gather_program_data(config.QUALTRICS_CSV_PATH)