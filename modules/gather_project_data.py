"""
gather_project_data.py
2024-01-03 ZD

This script defines primary function gather_project_data, which aggregates 
grants data by project_id and pulls relevant data from the grants to build a 
projects dataframe. 

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
    inconsistent_rows = df_lower[df_lower.groupby('queried_project_id')[field_name].transform('nunique') > 1].copy()
    inconsistent_rows['mismatched_field'] = field_name
    return inconsistent_rows

def get_newest_or_oldest_value(df, field_name, newest=True):
    """Returns the newest or oldest value of a field based on 'award_notice_date'."""

    return df.sort_values('award_notice_date', ascending=not newest).groupby('queried_project_id')[field_name].first()

def gather_project_data(grants_df):
    """Aggregates grants data into project-level information, handling inconsistencies and validation."""

    field_aggregation_map = {
        'queried_project_id': 'Identical',
        'project_title': 'Newest',
        'abstract_text': 'Newest',
        'project_start_date': 'Oldest',
        'project_end_date': 'Newest',
        'opportunity_number': 'Identical',
        'api_source_search': 'Identical',
        'org_name': 'Identical',
        'org_city': 'Identical',
        'org_state': 'Identical',
        'org_country': 'Identical',
        'program.program_id': 'Identical'
    }

    projects_df = pd.DataFrame(columns=field_aggregation_map.keys())

  # Validate and aggregate "Identical" fields, collecting mismatched rows
    df_mismatch = pd.DataFrame()
    for field_name, aggregation_type in field_aggregation_map.items():
        if aggregation_type == 'Identical':
            mismatched_rows = validate_identical_values(grants_df, field_name)
            df_mismatch = pd.concat([df_mismatch, mismatched_rows], ignore_index=True)
            projects_df[field_name] = grants_df.groupby('queried_project_id')[field_name].first()

    # Get "Newest" or "Oldest" fields
    for field_name, aggregation_type in field_aggregation_map.items():
        if aggregation_type in ('Newest', 'Oldest'):
            projects_df[field_name] = get_newest_or_oldest_value(grants_df, field_name, newest=aggregation_type == 'Newest')

    # Rename 'queried_project_id' to 'project_id'
    projects_df.rename(columns={'queried_project_id': 'project_id'}, inplace=True)

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