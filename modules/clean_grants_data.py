# clean_grants_data.py
# 2023-08-07 ZD

# This script defines primary function `clean_grants_data` intended to be run
# on grants data pulled from the NIH RePORTER API. It will process the grants 
# data into a format best suited for the Index of NCI Studies data ingestion
# process. 
#
# Cleaning Steps in order:
# 1. Remove unnecessary columns
# 2. Extract and format full names for PI and PO columns
# 3. Extract total NCI cost
# 4. Rename columns to match existing INS terms


import pandas as pd
import config

def clean_grants_data(grants_data):
    """Create clean dataframes from NIH RePORTER API response JSONs."""

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

    # Step 4: Rename columns to match INS terms
    rename_dict = config.API_FIELD_RENAMER
    df_renamed = df_cleaned_funding.copy()
    df_renamed.rename(columns=rename_dict)

    return df_renamed


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

