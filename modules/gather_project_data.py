"""
gather_project_data.py
2024-01-03 ZD

This script defines primary function gather_project_data, which aggregates 
grants data by project_id and pulls relevant data from the grants to build a 
projects dataframe. Some project_ids will be repeated if they have different 
program_id values.

The output `projects_df` and exported project.csv contain columns:
    - project_title
    - abstract_text
    - project_start_date
    - project_end_date
    - opportunity_number
    - api_source_search
    - org_name
    - org_city
    - org_state
    - org_country
    - program.program_id
"""

import os
import sys
from datetime import datetime

import pandas as pd

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


def validate_identical_values(df, field_name):
    """Identifies rows with mismatches in 'field_name', ignoring case differences."""

    # Skip this step for api_source_search and program values
    skip_cols = ['api_source_search', 'program.program_id']
    if field_name in skip_cols:
        return None

    # Convert values to lowercase for comparison
    df_lower = df.copy()
    df_lower[field_name] = df_lower[field_name].str.lower()

    # Find rows with more than one value within a field for a given project
    inconsistent_rows = df_lower[df_lower.groupby(
                                    ['queried_project_id','program.program_id'])
                                    [field_name].transform('nunique') > 1].copy()
    
    # Gather mismatched rows for reporting
    inconsistent_rows['mismatched_field'] = field_name

    return inconsistent_rows



def get_newest_or_oldest_value(df: pd.DataFrame, field_name: str, agg_type: str):
    """
    Retrieves the newest or oldest value of a specified field within each 
    project group.

    Args:
        df (pd.DataFrame): The DataFrame containing the data.
        field_name (str): The name of the field to extract the value from.
        agg_type (str): Specifies whether to retrieve the newest or oldest value.
            Valid options are 'newest' and 'oldest'.

    Returns:
        str: The extracted newest or oldest value of the specified field for 
            each project.

    Raises:
        ValueError: If an invalid agg_type is provided.
    """

    if agg_type == 'oldest':
        get_oldest = True
    elif agg_type == 'newest':
        get_oldest = False
    else:
        raise ValueError(f"Invalid agg_type input: {agg_type}. Check Config.")

    # Determine the appropriate sorting column based on the field name
    if field_name in ('project_start_date', 'project_end_date'):
        sort_col = field_name
    else:
        # Default sorting column
        sort_col = 'award_notice_date'

    # Use fiscal_year as a fallback if the default sorting column contains NaNs
    if df.groupby('queried_project_id')[sort_col].apply(pd.isna).any():
        sort_col = 'fiscal_year'

    # Sort values, group by project, and extract the first value
    value = (
        df.sort_values(sort_col, ascending=get_oldest)
        .groupby('queried_project_id')[field_name]
        .first()
    )

    return value



def get_list_values(df, field_name):
    """Collects unique values within a field and lists them as a 
    semicolon-separated string."""

    # Filter for string values if only strings are valid
    if not pd.api.types.is_string_dtype(df[field_name]):
        df = df[df[field_name].apply(lambda x: isinstance(x, str))]

    # Handle missing values and force values to string before sorting
    value = df.groupby(['queried_project_id'])[field_name].unique().apply(
        lambda x: ';'.join([str(v) for v in sorted(x)]))

    return value



def gather_project_data(grants_df):
    """Aggregates grants data into project-level information, handling 
    inconsistencies and validation."""

    # Define how values are selected from grants to populate project values
    # program.program_id is excluded and added later
    field_agg_map = {
        'queried_project_id': 'identical',
        'project_title': 'newest',
        'abstract_text': 'newest',
        'project_start_date': 'oldest',
        'project_end_date': 'newest',
        'opportunity_number': 'list',
        'api_source_search': 'identical',
        'org_name': 'identical',
        'org_city': 'identical',
        'org_state': 'identical',
        'org_country': 'identical',
    }

    # Build empty projects df and mismatch df
    projects_df = pd.DataFrame(columns=field_agg_map.keys())
    df_mismatch = pd.DataFrame()

    # Validate and aggregate fields based on aggregation type
    for field_name, agg_type in field_agg_map.items():

        # Validate identical values
        if agg_type == 'identical':
            # Gather any mismatched values that should be identical
            mismatched_rows = validate_identical_values(grants_df, field_name)
            df_mismatch = pd.concat([df_mismatch, mismatched_rows], 
                                    ignore_index=True)
            # Default to using values from newest projects 
            projects_df[field_name] = get_newest_or_oldest_value(grants_df, 
                                                                 field_name, 
                                                                 'newest')

        # Use the newest or oldest award date to select the project value
        elif agg_type in ('newest', 'oldest'):
            projects_df[field_name] = get_newest_or_oldest_value(grants_df, 
                                                                 field_name, 
                                                                 agg_type)

        # List all distinct values within a project as the project value(s)
        elif agg_type == 'list':
            projects_df[field_name] = get_list_values(grants_df, field_name)

    # Reset index to prep for merging programs
    projects_df = projects_df.reset_index(drop=True)

    # Get new df of all unique project-program combinations
    project_program_df = grants_df[['queried_project_id']
                                   ].drop_duplicates().reset_index(drop=True)

    # Add program IDs so each project has a single row for each program
    projects_df = projects_df.merge(project_program_df,
                                    on='queried_project_id',
                                    how='outer')

    # Rename 'queried_project_id' to 'project_id'
    projects_df.rename(columns={'queried_project_id': 'project_id'}, 
                       inplace=True)

    # Export as csv
    projects_df.to_csv(config.PROJECTS_INTERMED_PATH, index=False)
    print(f"Success! Project data exported to {config.PROJECTS_INTERMED_PATH}")

    # Export DataFrame with mismatched data for manual review
    if not df_mismatch.empty:
        mismatched_path = config.MISMATCHED_PROJECT_VALUES_REPORT
        df_mismatch.to_csv(mismatched_path, index=False)
        print(f"DataFrame with mismatched rows exported to {mismatched_path}")

    return projects_df



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Load Grants data
    grant_filepath = config.GRANTS_INTERMED_PATH
    grants_df = pd.read_csv(grant_filepath)
    print(f"Grants data loaded from {grant_filepath}.")

    # Build Projects data from grants
    gather_project_data(grants_df)