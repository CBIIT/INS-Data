"""
package_output_data.py
2023-12-21 ZD

This script defines primary function package_output_data that applies final 
formatting, filtering, etc. to the intermediate data gathering outputs in order
to prepare them for INS ingestion. These outputs are saved as TSVs. 
"""

import os
import sys
import unicodedata
import hashlib

import pandas as pd
import re


# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


def add_type_column(df, datatype):
    """Add a column for datatype and fill with specified datatype."""

    df_with_type = df.copy()
    df_with_type['type'] = datatype

    return df_with_type



def reorder_columns(df, column_configs, datatype):
    """Reorder and validate columns of df using dict from config.py. 
    Report any dropped columns and throw error if defined column not found.
    """
    # Check if the datatype exists in column_configs
    if datatype not in column_configs:
        raise ValueError(f"Invalid datatype: {datatype}")

    # Extract relevant configuration for the datatype
    config = column_configs[datatype]

    # Extract relevant information from the configuration
    keep_and_rename = config.get('keep_and_rename', {})

    # Check if all desired columns are present before renaming
    if not set(keep_and_rename.keys()).issubset(set(df.columns)):
        missing_columns = set(keep_and_rename.keys()) - set(df.columns)
        raise ValueError(f"Fields missing from intermediate {datatype} data: "
                         f"{', '.join(missing_columns)}")

    # Rename columns using keep_and_rename
    df.rename(columns=keep_and_rename, inplace=True)

    # Define list of columns to keep
    columns_to_keep = list(keep_and_rename.values())

    # Drop any columns not part of keep_and_rename
    dropped_columns = set(df.columns) - set(columns_to_keep)
    if dropped_columns:
        print(f"Columns dropped from {datatype} output: "
              f"{', '.join(dropped_columns)}")

    # Keep and reorder columns
    df_reordered = df[columns_to_keep]

    return df_reordered



def validate_first_columns(df, column_configs, datatype):
    """Validate that type, node_id, and link_id (opt) are first columns."""

    # Check if the datatype exists in column_configs and then get it
    if datatype not in column_configs:
        raise ValueError(f"Invalid datatype: {datatype}")
    config = column_configs.get(datatype)

    # Get values from column configuration
    node_id = config.get('node_id')
    link_id = config.get('link_id', None)

    # If link_id exists, add it to the list of expected columns 
    expected_cols = ['type', node_id]
    if link_id:
        expected_cols.append(link_id)

    # Check the order of columns in df
    if list(df.columns[:len(expected_cols)]) != expected_cols:
        raise ValueError(f"First columns in {datatype} output are not in the "
                         f"expected order. Reorder columns in config.py to fix.\n"
                         f"Expected order: {', '.join(expected_cols)}.\n"
                         f"Actual order:   "
                         f"{', '.join(df.columns[:len(expected_cols)])}")
    
    return None



def replace_defined_characters(text):
    """Use a mappping to replace specific non-standard characters in text. Any
    non-ascii characters not included here will be dropped."""

    # Return NaN and non-string values as-is
    if pd.isna(text) or not isinstance(text, str):
        return text

    # Define translation table to replace characters
    translation_table = str.maketrans({
        '\"'    : '\'',     # STANDARD double quotes - replace with single quote
        '\u2028': ' ',      # Line Separator (LS) - replace with space
        '\u2029': ' ',      # Paragraph Separator (PS) - replace with space
        '\u00a0': ' ',      # Non-breaking space - replace with space
        '\ufeff': '',       # Byte Order Mark (BOM) warning sign - remove
        '\u2013': '-',      # En Dash - replace with dash
        '\u2014': '-',      # Em Dash - replace with dash
        '\u2010': '-',      # Hyphen - replace with dash
        '\u201c': '\'',     # Left double quotes - replace with single quote
        '\u201d': '\'',     # Right double quotes - replace with single quote
        '\u2018': '\'',     # Left single quote - replace with single quote
        '\u2019': '\'',     # Right single quote - replace with single quote
        '\u0004': '',       # Unused unicode - remove
        '\u0005': '',       # Unused unicode - remove
        '\u0006': '',       # Unused unicode - remove
        '\u0007': '',       # Unused unicode - remove
        '\u03b1': 'a',      # Greek alpha - replace with a
        '\u03b2': 'B',      # Greek beta - replace with B
        '\u03b3': 'g',      # Greek gamma - replace with g
        '\u03b4': 'd',      # Greek lower delta
        '\u03F5': 'e',      # Greek lower epsilon
        '\u03ba': 'k',      # Greek kappa - replace with k
        '\u03bc': 'u',      # Greek mu - replace with u
        '\u03BB': 'l',      # Greek lambda - replace with l
        '\u0394': 'D',      # Greek upper delta
        '\u2081': '1',      # Subscript 1
        '\u2082': '2',      # Subscript 2
        '\u2083': '3',      # Subscript 3
        '\u2084': '4',      # Subscript 4
        '\u2085': '5',      # Subscript 5
        '\u2086': '6',      # Subscript 6
        '\u2087': '7',      # Subscript 7
        '\u2088': '8',      # Subscript 8
        '\u2089': '9',      # Subscript 9
        '\u00a9': '(C)',    # Copywrite symbol
        '\u2264': '<=',     # Less than or equal to
        '\u2265': '>=',     # Greater than or equal to

        # Add mappings here in the future as needed
        })

    # Use translate only for strings
    replaced_text = (text.translate(translation_table) 
                     if isinstance(text, str) else text)
    
    # Use translate only for strings
    return replaced_text



def normalize_encoding(text):
    """Normalizes text encoding for consistency and compatibility."""

    # Return NaN and non-string values as-is
    if pd.isna(text) or not isinstance(text, str):
        return text

    # Encode and decode text through an ascii "filter" to normalize
    try:
        # Using NFD will decompose characters into base and combining accent
        formatted_text = (unicodedata.normalize('NFD', text)
                # 'ignore' argument will drop any non-ascii, including accents
                          .encode('ascii', 'ignore').decode('utf-8'))
        return formatted_text
    
    # Second catch to handle non-string values as-is
    except TypeError:
        return text



def process_special_characters(df):
    """Cleans a DataFrame by normalizing special characters and encoding."""

    # Replace specific non-standard characters
    df = df.map(replace_defined_characters)
    # General normalization to utf-8 encoding
    df_cleaned = df.map(normalize_encoding)

    return df_cleaned



def remove_html_tags(text):
    """Removes HTML tags using regular expressions."""

    # Matches opening or closing tags with tag name like <i>italics</i>
    tag_pattern = r'<(/?)(\w+)>'

    return re.sub(tag_pattern, '', text)



def remove_html_tags_from_df(df, column_configs, datatype):
    """Removes HTML tags such as <i>italics</i> or <b>bold</b> from text."""

    # Get predefined columns containing HTML tags from configuration
    html_tag_cols = column_configs[datatype].get('html_tag_cols', None)

    # Iterate through expected HTML-tagged columns
    if html_tag_cols is not None:

        for col in html_tag_cols:
        
            # Check that expected column exists
            if col not in df.columns:
                raise ValueError(f"Expected html-tagged column {col} not found "
                                 f"within data. Check data or redefine config.")
        
            # Use regex to remove HTML tags but keep enclosed text
            df[col] = df[col].apply(remove_html_tags)

    return df



def format_datetime_columns(df, column_configs, datatype):
    """Converts datetime columns to "yyyy-mm-dd" strings."""

    # Get predefined list-like columns from configuration
    datetime_cols = column_configs[datatype].get('datetime_cols', None)

    # Iterate through expected datetime-like columns
    if datetime_cols is not None:

        for col in datetime_cols:

            # Check that expected column exists
            if col not in df.columns:
                raise ValueError(f"Expected datetime column {col} not found "
                                 f"within data. Check data or redefine config.")
            
            # Convert datetime string values to datetime and back to string
            format_string = '%Y-%m-%d'
            df[col] = (pd.to_datetime(df[col], utc=True)
                       .dt.strftime(format_string))

    return df



def validate_listlike_columns(df, column_configs, datatype):
    """Validate semicolon separation without spaces in list-like columns."""

    # Get predefined list-like columns from configuration
    listlike_cols = column_configs[datatype].get('list_like_cols', None)

    # Iterate through expected list-like columns
    if listlike_cols is not None:

        for col in listlike_cols:

            # Check that expected column exists
            if col not in df.columns:
                raise ValueError(f"Expected list-like column {col} not found "
                                 f"within data. Check data or redefine config.")

            # Replace any semicolons followed by whitespace with just semicolon
            df[col] = df[col].str.replace('; ', ';')

            # Check if any semicolons exist in any value of the column
            if not df[col].astype(str).str.contains(';').any():
                print(f"WARNING: No semicolons found in list-like column {col}")

    return df



def validate_and_clean_unique_nodes(df, column_configs, datatype):
    """
    Validate that all node_ids are either unique values or valid duplicates.

    Args:
        df (pd.DataFrame): Input DataFrame
        column_configs (dict): Configuration settings for columns, including 
            node_id and link_id
        datatype (str): Type of data in the DataFrame. Indicates which info
            will be accessed from column_configs

    Returns:
        pd.DataFrame: A new DataFrame with invalid duplicate rows removed.

    Raises:
        ValueError: If the specified node_id column is not found in the DataFrame.

    Notes:
        This function checks for valid uniqueness or duplicates in node_ids 
            based on the provided configuration.

        - For many-to-many nodes, it ensures that any duplicates only differ 
            in the link_id column.
        - True duplicates (identical values in all columns) are kept in the 
            first row for production.
    """

    # Get node_id and link_id from column configurations
    node_id = column_configs[datatype].get('node_id')
    link_id = column_configs[datatype].get('link_id', None)

    # Check if the specified node_id column exists
    if node_id not in df.columns:
        raise ValueError(f"Node ID column '{node_id}' not found in DataFrame.")

    # Get only rows with duplicated node_ids 
    dup_nodes = df[df.duplicated(subset=node_id, keep=False)]
    
    # Get true duplicates where values in all columns are identical
    true_duplicates = dup_nodes[dup_nodes.duplicated(keep=False)].assign(
        error_reason=f"True duplicate in all columns. "
                     f"Kept first row for production.")
    
    # Empty dataframe to store mismatch rows
    mismatch_df = pd.DataFrame()

    # Get list of value columns to check for mismatch
    value_cols = list(set(df.columns) - set([node_id, link_id]))

    # Group rows by node_id for shared value checking
    for node, node_df in dup_nodes.groupby(node_id):
        
        # Check for mismatched values in other columns within the same node
        value_mismatches = node_df[~node_df.duplicated(
            subset=value_cols, keep=False)
            ].assign(error_reason=f"Mismatched values in columns with same "
                                  f"node_id. Both dropped for production.")

        # Add to running list of rows with mismatched values
        mismatch_df = pd.concat([mismatch_df, value_mismatches])

    # Combine list of annotated invalid duplicates for export to csv
    invalid_duplicates_df = pd.concat([true_duplicates, mismatch_df])

    # Save any invalid duplicate rows to a report csv for troubleshooting
    if len(invalid_duplicates_df) > 0:
        invalid_duplicates_path = f"{config.REMOVED_DUPLICATES}{datatype}.csv"
        os.makedirs(os.path.dirname(invalid_duplicates_path), exist_ok=True)
        invalid_duplicates_df.to_csv(invalid_duplicates_path, index=False)
        print(f"Invalid duplicate rows cleaned from output and recorded in "
              f"{invalid_duplicates_path}.")

    # Return a cleaned df with only unique node_ids or valid duplicates
    # Remove all mismatched node_ids and any true duplicates after the first
    if not mismatch_df.empty:
        valid_df = df[~df[node_id].isin(mismatch_df[node_id])
                    & ~df.duplicated(keep='first')
                    ].reset_index(drop=True)
        
    # If no mismatches found, then just remove true duplicates
    else:
        valid_df = df[~df.duplicated(keep='first')
                       ].reset_index(drop=True)

    return valid_df



def standardize_data(df, column_configs, datatype):
    """Group standardization functions common to all data types"""

    # Edit data to standardize
    df = add_type_column(df, datatype)
    df = reorder_columns(df, column_configs, datatype)
    df = process_special_characters(df)
    df = remove_html_tags_from_df(df, column_configs, datatype)
    df = validate_listlike_columns(df, column_configs, datatype)
    df = format_datetime_columns(df, column_configs, datatype)
    df = validate_and_clean_unique_nodes(df, column_configs, datatype)

    # Validate that data meets loading standards
    validate_first_columns(df, column_configs, datatype)
    
    return df



def package_programs(df_programs, column_configs):
    """Package programs data for INS loading."""

    print(f"---\nFinalizing TSV for program data...")

    # Standardize and validate data
    df_programs_output = standardize_data(df_programs, column_configs, 
                                        datatype='program')

    # Export as TSV
    output_filepath = config.PROGRAMS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_programs_output.to_csv(output_filepath, sep='\t', index=False, 
                            encoding='utf-8')

    print(f"Done! Final program data saved as {output_filepath}.")

    return df_programs_output



def package_grants(df_grants, column_configs):
    """Package grants data for INS loading."""

    print(f"---\nFinalizing TSV for grant data...")

    # Standardize and validate data
    df_grants_output = standardize_data(df_grants, column_configs, 
                                        datatype='grant')

    # Export as TSV
    output_filepath = config.GRANTS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_grants_output.to_csv(output_filepath, sep='\t', index=False, 
                            encoding='utf-8')

    print(f"Done! Final grant data saved as {output_filepath}.")

    return df_grants_output



def package_projects(df_projects, column_configs):
    """Package projects data for INS loading."""

    print(f"---\nFinalizing TSV for project data...") 

    # Standardize and validate data
    df_projects_output = standardize_data(df_projects, column_configs, 
                                        datatype='project')

    # Export as TSV
    output_filepath = config.PROJECTS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_projects_output.to_csv(output_filepath, sep='\t', index=False, 
                            encoding='utf-8')

    print(f"Done! Final projects data saved as {output_filepath}.")

    return df_projects_output



def package_publications(df_publications, column_configs):
    """Package publications data for INS loading."""

    print(f"---\nFinalizing TSV for publication data...") 

    # Standardize and validate data
    df_publications_output = standardize_data(df_publications, column_configs, 
                                        datatype='publication')

    # Export as TSV
    output_filepath = config.PUBLICATIONS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_publications_output.to_csv(output_filepath, sep='\t', index=False, 
                                    encoding='utf-8')

    print(f"Done! Final publication data saved as {output_filepath}.")

    return df_publications_output



def remove_publications_before_projects(df_publications: pd.DataFrame, 
                                        df_projects: pd.DataFrame, 
                                        day_diff_allowed:int=0) -> pd.DataFrame:
    """Remove publications published before the associated project start date.
    
    Args:
        df_publications (pd.DataFrame): DataFrame containing publication data.
        df_projects (pd.DataFrame): DataFrame containing project information.
        day_diff_allowed (int, optional): Allowable difference in days between
            publication date and project start date. Defaults to 0. 
            If None, then no filtering is applied and original publications
            are returned.

    Returns:
        pd.DataFrame: DataFrame containing publications with valid dates. 
            If day_diff_allowed == None, returns original df_publications.

    Raises:
        ValueError: If column structure of the filtered DataFrame is altered

    Notes:
        Exports the DataFrame of removed publications as CSV to a location
            defined in config.py.
    """

    # Skip this function entirely if None provided for day_diff_allowed 
    if not day_diff_allowed:
        print(f"Filter NOT applied for publications published before project "
              f"start date. Keeping all. Change this in config.py if desired.")
        return df_publications

    # Validate day difference and convert to pd.Timedelta value
    if not isinstance(day_diff_allowed, int):
        raise ValueError(f"Expected int value for `day_diff_allowed` but found "
                         f"{type(day_diff_allowed)}. Fix this in config.py.")
    diff = str(day_diff_allowed) + ' days'

    # Merge project start dates with publication data
    df_merged = pd.merge(left=df_publications, 
                        right=df_projects[['project_id','project_start_date']], 
                        how='left',left_on='coreproject', right_on='project_id',
                        # Ignore duplicate projects listed for multiple programs
                        ).drop_duplicates()

    # Get subset where publication pub date is earlier than project start date
    df_early_pubs = df_merged[
        # Convert dates to datetime for more accurate comparison
        (pd.to_datetime(df_merged['publication_date'],utc=True) 
         # Add in allowable day difference
        + pd.Timedelta(diff))
        < pd.to_datetime(df_merged['project_start_date'],utc=True)
        # Set on a copy to avoid SettingWithCopyWarning
        ].copy()

    # Get publications to keep after early pubs removed
    df_keep_pubs = df_merged[~df_merged.index.isin(df_early_pubs.index)]

    # Drop merged columns and reindex filtered publications
    df_keep_pubs = df_keep_pubs.drop(columns=['project_id','project_start_date']
                                     ).reset_index(drop=True)
    
    # Add column showing difference between pub date and project start date
    df_early_pubs.loc[:,'day_diff'] = (
        pd.to_datetime(df_early_pubs['project_start_date'],utc=True) -
        pd.to_datetime(df_early_pubs['publication_date'],utc=True)
        ).dt.days

    
    if not df_early_pubs.empty:

        # Define and make directory for report
        early_pub_report_path = config.REMOVED_EARLY_PUBLICATIONS
        os.makedirs(os.path.dirname(early_pub_report_path), exist_ok=True)

        # Export report of removed early publications
        df_early_pubs.to_csv(early_pub_report_path, index=False)
        print(f"---\nEarly Publications detected:\n"
              f"{len(df_early_pubs)} Publications with publication date more "
              f"than {diff} before the associated project start date were "
              f"removed and saved to {early_pub_report_path}")
        
    # Validate that columns have not changed in returned filtered list
    if df_publications.columns.tolist() != df_keep_pubs.columns.tolist():
        raise ValueError(f"Unexpected change in columns during "
                         f"`remove_publications_before_projects()`.\n"
                         f"Expected: {df_publications.columns.tolist()}\n"
                         f"Actual:   {df_keep_pubs.columns.tolist()}")

    return df_keep_pubs



def get_enum_values(df, cols, output="enums.txt"):
    """Extracts unique values from specified columns in a DataFrame,
    sorts them, and writes them to a YAML-like text file.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        cols (list): A list of column names to extract unique values from.
        output (str, optional): The output file path. Defaults to "enums.txt".
    """

    # Make directory if it doesn't already exist
    os.makedirs(os.path.dirname(output), exist_ok=True)

    # Get data type
    df_type = df['type'][0]

    # Build empty dict
    enums_dict = {}

    for col in cols:
        # Get unique values from semicolon-separated strings
        unique_vals = df[col].dropna().str.split(";").explode().unique()
        # Sort for conistency
        unique_vals = sorted(list(unique_vals))
        # Add "Enum" header to match data model yaml structure
        enums_dict[col] = {"Enum": unique_vals}

    # Save as txt with yaml-like formatting 
    with open(output, "w") as f:
        # Use 2-space indents
        indent = "  "

        # Manually imitate formatting of yaml structure and get data type
        f.write(f"PropDefinitions:\n"
                f"{indent}#properties of {df_type}\n"
                f"{indent}...\n")
        
        # Add columns as headers
        for col, enum_data in enums_dict.items():
            f.write(f"{indent}{col}:\n"
                        f"{indent*2}... \n"
                        f"{indent*2}Type:\n"
                                f"{indent*3}value_type: list\n"
                                f"{indent*3}Enum:\n")
            
            for val in enum_data["Enum"]:
                            # Iterate through unique values as enums
                            f.write(f"{indent*4}- {val}\n")
            # Add trailing ... to end of section
            f.write(f"{indent*2}...\n")

    print(f"Done! Enumerated values for {df_type}s: {cols} saved to {output}")



def generate_md5_hash(file_path):
    """Generates the MD5 hash for a given file."""

    with open(file_path, "rb") as f:
        data = f.read()
        md5_hash = hashlib.md5(data).hexdigest()

    return md5_hash



def generate_md5_hash_file(directory):
    """Generate MD5s for all TSVs in a directory and list in a text file."""

    hashes_file = os.path.join(directory, "_md5.txt")

    with open(hashes_file, "w") as f:
        for root, _, files in os.walk(directory):
            for file in files:
                if file.endswith(".tsv"):
                    file_path = os.path.join(root, file)
                    md5_hash = generate_md5_hash(file_path)
                    f.write(f"{md5_hash}\t{file_path}\n")
    
    print(f"Done! MD5 hashes saved to {hashes_file}.")



def package_output_data():
    """Run all data packaging steps for all data types."""

    print(f"\n---\nDATA PACKAGING:\n"
          f"Performing final data packaging steps...\n---\n")

    # Single data model dict defined in config
    column_configs = config.COLUMN_CONFIGS 

    # Load programs data
    if os.path.exists(config.PROGRAMS_INTERMED_PATH):
        programs_exist = True
        df_programs = pd.read_csv(config.PROGRAMS_INTERMED_PATH)
        print(f"Loaded Programs file from {config.PROGRAMS_INTERMED_PATH}")
    else: 
        programs_exist = False
        print(f"No Programs file found.")
    
    # Load grants data
    if os.path.exists(config.GRANTS_INTERMED_PATH):
        grants_exist = True
        df_grants = pd.read_csv(config.GRANTS_INTERMED_PATH)
        print(f"Loaded Grants file from {config.GRANTS_INTERMED_PATH}")
    else: 
        grants_exist = False
        print(f"No Grants file found.")

    # Load projects data
    if os.path.exists(config.PROJECTS_INTERMED_PATH):
        projects_exist = True
        df_projects = pd.read_csv(config.PROJECTS_INTERMED_PATH)
        print(f"Loaded Projects file from {config.PROJECTS_INTERMED_PATH}")
    else: 
        projects_exist = False
        print(f"No Projects file found.")

    # Load publications data
    if os.path.exists(config.PUBLICATIONS_INTERMED_PATH):
        publications_exist = True
        df_publications = pd.read_csv(config.PUBLICATIONS_INTERMED_PATH)
        print(f"Loaded Publications file from {config.PUBLICATIONS_INTERMED_PATH}")
    else: 
        publications_exist = False
        print(f"No Publications file found.")

    # Special handling
    print(f"---\nApplying special handling steps...")
    if publications_exist and projects_exist:
        df_publications = remove_publications_before_projects(
                                df_publications,
                                df_projects,
                                day_diff_allowed=config.PUB_PROJECT_DAY_DIFF)

    # Final packaging
    if programs_exist:
        df_programs_out = package_programs(df_programs, column_configs)
    if grants_exist:
        df_grants_out = package_grants(df_grants, column_configs)
    if projects_exist:
        df_projects_out = package_projects(df_projects, column_configs)
    if publications_exist:
        df_publications_out = package_publications(df_publications, column_configs)

    # Special post-processing handling
    print(f"---\nGenerating enumerated values for data model...")
    if programs_exist:
        get_enum_values(df_programs_out, 
                        cols = config.ENUM_PROGRAM_COLS, 
                        output = config.ENUM_PROGRAM_PATH)
        
    # Generate md5.txt
    print(f"---\nGenerating md5 hashes for file validation...")
    generate_md5_hash_file(config.OUTPUT_GATHERED_DIR)

# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Run main packaging function for all data
    package_output_data()