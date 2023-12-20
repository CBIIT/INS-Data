"""
gather_grants_data.py
2023-08-07 ZD

This script defines primary function gather_grants_data. 

First, get_nih_reporter_grants calls the public NIH RePORTER API to gather 
grants data for the provided NOFO (Notice of Funding Opportunity) and/or awards 
(i.e. grants, supplements, or parent projects) for each Key Program in the Key 
Programs CSV. Next, clean_grants_data calls helper functions to process the 
grants data to filter and reformat for use within INS.

NIH RePORTER Endpoint: https://api.reporter.nih.gov/v2/projects/search 
NIH RePORTER Docs: https://api.reporter.nih.gov/

Cleaning Steps in order:
1. Remove unnecessary columns
2. Extract and format full names for PI and PO columns
3. Extract total NCI cost and remove grants without NCI funding
4. Extract organization information
5. Rename columns to match existing INS terms
6. Drop duplicate grant rows
"""

import os
import sys
from datetime import datetime
from time import sleep  # for retrying API calls
from math import ceil  # for pagination logging

import pandas as pd
import re
import requests
from tqdm import tqdm # for progress bars

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


def get_nih_reporter_grants(search_values:str, 
                            search_type:str, 
                            print_meta=False):
    """Get grants info from NIH RePOTER for either NOFOs or Awards as specified.

    Args:
        search_values (str): String of values to search, separated by semicolons
        search_type (str): The type of search, e.g., 'nofo' or 'award'
        print_meta (bool, optional): A boolean indicator. If True, print API 
        gathering process results to the console. Defaults to False.

    Returns:
        grants_data (str): JSON API response containing grants data
    """

    # Set the search field based on the search_type
    if search_type == 'nofo':
        search_field = "opportunity_numbers"
    elif search_type == 'award':
        search_field = "project_nums"
    else:
        raise ValueError("Invalid search type.")

    base_url = "https://api.reporter.nih.gov/v2/projects/search"
    grants_data = []

    for search_value in search_values:
        # Check for blank search values
        if not search_value:
            if print_meta == True:
                tqdm.write(f"Blank {search_type} value encountered. Skipping.")
            continue

        # Set default values for params not likley to change
        LIMIT = 500
        MAX_ATTEMPTS = 5
        RETRY_TIME = 2
        EARLIEST_FISCAL_YEAR = config.API_EARLIEST_FISCAL_YEAR

        # Get list of fiscal years to query
        current_year = datetime.now().year
        fiscal_year_list = [str(year) for year in range(
            EARLIEST_FISCAL_YEAR, current_year + 1)]

        # Set starting value for counters
        offset = 0
        page = 0
        attempts = 0

        # RePORTER API sets a max limit of 500 records per call.
        # Keep looping each call in "pages" until the number of records 
        # gathered reaches the total number of records available. 
        
        # Set a cap on the number of attempts at a failed call before 
        # moving on to the next award.
        while attempts < MAX_ATTEMPTS:
            # Set parameters for API call
            params = {
                "criteria": {
                    search_field: [search_value],
                    "exclude_subprojects": True,
                    "agencies": ["NCI"],
                    "is_agency_funding": True,
                    "fiscal_years": fiscal_year_list
                },
                "limit": LIMIT,
                "offset": offset,
                "sort_field": "appl_id",
                "sort_order": "desc",
            }

            try: 
                # Define response details
                response = requests.post(base_url, 
                                        json=params, 
                                        headers={
                                            "accept": "application/json", 
                                            "Content-Type": "application/json"})

                # If response is good, get results
                if response.status_code == 200:
                    grants = response.json()
                    # Add API source indicator
                    for grant in grants['results']:
                        grant['api_source_search'] = f"{search_type}_{search_value}"
                    # Add grants to running list
                    grants_data.extend(grants['results'])

                    # Increase offset by limit to get next "page"
                    total_records = grants['meta']['total']
                    offset = offset + LIMIT
                    page = page + 1

                    # Print paginated partial optional metadata
                    # Consider replacing this with proper logging
                    if print_meta == True:
                        total_pages = max(ceil(total_records/LIMIT),1)
                        tqdm.write(f"{search_type}: {search_value} "
                              f"({page}/{total_pages}): {grants['meta']}")

                    # Stop looping if offset has reached total record count
                    if offset >= total_records:
                        break

                # Handle 500 errors by retrying after 2 second delay
                elif response.status_code == 500:
                    attempts = attempts + 1
                    tqdm.write(f"Received a 500 error for "
                          f"{search_type} '{search_value}'. "
                          f"Retrying after {RETRY_TIME} seconds. "
                          f"Attempt {attempts}/{MAX_ATTEMPTS}")
                    sleep(RETRY_TIME)
                else:
                    tqdm.write(f"Error occurred while fetching grants for "
                          f"{search_type} '{search_value}': "
                          f"{response.status_code}")
                    break

            except requests.exceptions.RequestException as e:
                tqdm.write(f"An error occurred while making the API call for "
                      f"{search_type} '{search_value}': {e}")
                break

    return grants_data



def concatenate_full_names(row):
    """Replace JSON row values with full names list from within JSON."""

    # Check for blank row, then get each full name
    if row:
        full_names = [item['full_name'] for item in row]
        # Concatenate with a comma space between
        return ', '.join(full_names)
    else:
        # Return blank output for blank input
        return ''
    
    

def format_name_column(name_str):
    """Format capitalization and spacing of full names."""

    # Capitalize the first letter of each name
    formatted_name = name_str.title()
    # Remove double whitespaces between names
    formatted_name = formatted_name.replace('  ',' ')
    return formatted_name



def extract_total_cost(fundings):
    """Extract the total NCI funding from IC award totals"""

    for funding in fundings:
        # CA is the NCI code
        if funding['code'] == 'CA':
            return int(funding['total_cost'])
    # Return 0 if no NCI funding found
    return int(0)



def clean_abstract(text):
    """Clean abstract text of odd characters from mismatched encoding"""

    # Define a list of characters to be removed or replaced
    chars_to_remove = ['\xad']
    chars_to_space = ['  ']

    if isinstance(text, str):
        # Remove unwanted characters (no space added)
        for char in chars_to_remove:
            text = text.replace(char, '')
        # Replace unwanted characters with a space
        for char in chars_to_space:
            text = text.replace(char, ' ')

        # Remove non-breaking spaces, newlines, and other unwanted characters
        cleaned_text = re.sub(r'[\s\xa0]+', ' ', text).strip()
        return cleaned_text
        
    # Return original value if it is NaN or float
    else: 
        return text



def format_organization_columns(df, org_field_old, org_subfields):
    """Extract organization fields from nested API JSON organization field.

    Args:
        df (pd.DataFrame): pandas DataFrame containing columns with nested data
        org_field_old (str): The name of the old nested API JSON organization 
            field
        org_subfields (list): A list of strings representing the organization 
            subfields to extract

    Returns:
        pandas.DataFrame: The DataFrame with new columns for each extracted 
            organization subfield.
    """

    # Extract defined org subfields as new columns and fill with values
    for org_subfield in org_subfields:
        df[org_subfield] = (df[org_field_old]
                            .apply(lambda row: row[org_subfield]))
    # Drop the old nested organization column after extracting contents
    df.drop(org_field_old, axis=1, inplace=True)
    
    return df



def clean_grants_data(grants_data, print_meta=False):
    """Create clean DataFrames from NIH RePORTER API response JSONs.

    This function performs several steps to clean and process the input data:
    
    1. Load JSON grants data as a pandas DataFrame and select specified columns
    2. Extract and format Principal Investigator (PI) and Program Officer (PO) 
        full names from nested JSON
    3. Extract total National Cancer Institute (NCI) cost from nested JSON field
    4. Remove grants that did not receive NCI funding
    5. Extract desired organization values from nested JSON field
    6. Rename columns to match INS terms
    7. Clean abstract text
    8. Drop duplicate grant rows

    Args:
        grants_data (list): List of NIH RePORTER API response JSONs
        print_meta (bool): Boolean indicator. If True, print process results to 
            the console.

    Returns:
        pandas.DataFrame: Cleaned DataFrame containing grants information

    Notes:
        Exports a project.csv to the versioned intermediate directory
        This function uses various helper functions and configuration parameters
            defined in the config module
    """

    # Step 1: Load JSON grants data as pandas dataframe and select columns
    kept_fields = config.API_FIELDS
    df_grants = pd.DataFrame(grants_data)[kept_fields]

    # Step 2: Extract and format PI and PO full names from nested JSON
    name_cols = config.API_NAME_FIELDS
    df_cleaned_names = df_grants.copy()
    df_cleaned_names[name_cols] = (df_cleaned_names[name_cols]
                                   .apply(lambda col: col
                                          .apply(concatenate_full_names)
                                          .apply(format_name_column)))

    # Step 3: Extract total NCI cost from nested JSON field
    agency_funding_col = config.API_AGENCY_FUNDING_FIELD
    df_cleaned_funding = df_cleaned_names.copy()
    df_cleaned_funding[agency_funding_col] = (df_cleaned_funding
                                              [agency_funding_col]
                                              .apply(extract_total_cost))
    
    # Step 4: Remove grants that did not receive NCI funding
    df_cleaned_funding_nci_only = (df_cleaned_funding[df_cleaned_funding
                                    [agency_funding_col] > 0]
                                    .reset_index(drop=True))
    # Print output to show number of grants removed
    if print_meta == True:
        all_grant_count = len(df_cleaned_funding)
        nci_grant_count = len(df_cleaned_funding_nci_only)
        if all_grant_count > nci_grant_count:
            tqdm.write(f"{all_grant_count-nci_grant_count} grants without NCI "
                f"funding removed from list. \n"
                f"{nci_grant_count} NCI-funded grants remain.")
    
    # Step 5: Extract desired organization values from nested JSON field
    org_field_old = config.API_ORG_FIELD
    org_fields_keep = config.API_ORG_SUBFIELDS
    df_cleaned_orgs = df_cleaned_funding_nci_only.copy()
    df_cleaned_orgs = (format_organization_columns(df_cleaned_orgs,
                                                   org_field_old,
                                                   org_fields_keep))

    # Step 6: Rename columns to match INS terms
    rename_dict = config.API_FIELD_RENAMER
    df_renamed = df_cleaned_orgs.copy()
    df_renamed.rename(columns=rename_dict, inplace=True)

    # Step 7: Clean abstract text
    abstract_col = config.ABSTRACT_TEXT_FIELD
    df_cleaned_abstract = df_renamed.copy()
    df_cleaned_abstract[abstract_col] = (df_cleaned_abstract
                                         [abstract_col]
                                         .apply(clean_abstract))
    
    # Step 8: Drop duplicate grant rows
    df_duplicates_removed = df_cleaned_abstract.copy()
    df_duplicates_removed.drop_duplicates(inplace=True, ignore_index=True)
    # Print output to show number of duplicates removed
    if print_meta == True:
        if len(df_cleaned_abstract) > len(df_duplicates_removed):
            tqdm.write(f"{len(df_cleaned_abstract)-len(df_duplicates_removed)} "
                f"duplicates removed. \n"
                f"{len(df_duplicates_removed)} grants remain.")

    return df_duplicates_removed



def gather_grants_data(key_programs_df, print_meta=False):
    """Use the NIH RePORTER to gather grants data for each Program and then 
    process the grants data for use within INS.

    Args:
        key_programs_df (pd.DataFrame): DataFrame containing Program info. Must 
            include 'nofo' and/or 'award' columns listing API search terms. 
        print_meta (bool): Boolean indicator. If True, print API gathering 
            process results to the console.

    Returns:
        pd.DataFrame: Cleaned DataFrame containing grants information for each 
            Key Program.
    """

    print(f"---\nGathering, cleaning, and saving grants data from NIH "
        f"Reporter API for each Key Program...")

    # Create empty DataFrame to fill with grants
    all_cleaned_grants = pd.DataFrame()

    # Initialize a counter for total records
    total_records_count = 0

    # Iterate through each Key Program to get name, NOFOs, and Awards
    for index, row in tqdm(key_programs_df.iterrows(), 
                           total=len(key_programs_df), ncols=80):
        program_name = row['program_name']
        nofo_str = row['nofo']
        award_str = row['award']

        if print_meta == True:
            tqdm.write(f"---\n{program_name}")

        # Convert semicolon-separated award string into listed values
        if isinstance(award_str, str):
            award_list = [x.strip() for x in award_str.split(';')]
        else:
            award_list = []
        # Convert semicolon-separated NOFO string into listed values
        if isinstance(nofo_str, str):
            nofo_list = [x.strip() for x in nofo_str.split(';')]
        else:
            nofo_list = []

        # If neither NOFO nor Awards provided
        if pd.isna(nofo_str) and pd.isna(award_str):
            tqdm.write(f"NOTICE: No NOFOs or Awards defined for {program_name}." 
                       f"Skipping program.")
            continue

        # Gather grants data from NIH RePORTER
        award_grants_data = get_nih_reporter_grants(award_list, 'award', 
                                                    print_meta=print_meta)
        nofo_grants_data = get_nih_reporter_grants(nofo_list, 'nofo', 
                                                   print_meta=print_meta)

        # Combine and report
        grants_data = award_grants_data + nofo_grants_data
        total_records_count = total_records_count + len(grants_data)
        if print_meta == True:
            tqdm.write(f"{len(grants_data)} grants found.")
        if award_grants_data and nofo_grants_data:
            tqdm.write(f"NOTICE: Both Awards and NOFOs provided. Results may "
                       f"include duplicate grants gathered by both NOFO and "
                       f"Award searches.")
            continue

        # Skip remaining steps if no grants data gathered
        if len(grants_data) == 0:
            tqdm.write(f"NOTICE: No grants found for {program_name}")
            continue

        # Data Cleaning - Clean and format grants data
        cleaned_grants_data = clean_grants_data(grants_data, 
                                                print_meta=print_meta)

        # Add program id column from alphanumeric program name
        program_id = ''.join(filter(str.isalnum, program_name))
        cleaned_grants_data[config.PROGRAM_ID_FIELDNAME] = program_id

        # Combine and save cleaned grants data
        all_cleaned_grants = pd.concat([all_cleaned_grants,
                                        cleaned_grants_data]) 

    # Define versioned output directory using config.py
    project_filepath = config.PROJECTS_INTERMED_PATH
    os.makedirs(os.path.dirname(project_filepath), exist_ok=True)

    # Sort by program and project for consistency
    all_cleaned_grants.sort_values(by=[config.PROGRAM_ID_FIELDNAME,
                                    config.PROJECT_ID_FIELDNAME], 
                                    inplace=True, ignore_index=True)

    # Export to csv
    all_cleaned_grants.to_csv(project_filepath, index=False)
    
    print(f"---\nSuccess! NIH RePORTER API data gathered, cleaned, and saved.\n"
        f"{total_records_count} grants gathered across all Awards and NOFOs.\n"
        f"{len(all_cleaned_grants)} NCI-funded grants retained for INS. \n"
        f"Results can be found in {project_filepath}.\n---") 

    return all_cleaned_grants



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Load cleaned programs data
    key_programs_df = pd.read_csv(config.CLEANED_KEY_PROGRAMS_CSV)

    # Gather grants data
    all_cleaned_grants = gather_grants_data(key_programs_df, print_meta=False)