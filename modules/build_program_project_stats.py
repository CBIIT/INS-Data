"""
build_program_project_stats.py
2023-09-08 ZD

This script defines the main function build_program_project_stats that
generates reports summarizing program and project data gathered for INS.
These reports are not ingested into the site, but are useful for testing,
validation, and quick feedback during the pipeline.

Reports are output to the reports/ directory with a versioning structure
identical to the data/ directory.
"""

import pandas as pd
import os
from itertools import combinations

import config


def get_grant_stats_by_program(all_grants_data:pd.DataFrame):
    """Group grants by program and get aggregated summaries"""

    # Summarize columns of interest using aggregations defined in config
    grant_stats_by_program = (all_grants_data.groupby(
                                    config.PROGRAM_ID_FIELDNAME)
                                    .agg(config.STAT_AGG_FUNCS_BY_COL)
                                    .reset_index())
    
    # Rename fiscal year stat column for clarity
    fiscal_year_renamed = 'earliest_'+ config.STAT_FISCALYEAR_COL
    grant_stats_by_program.rename(columns={
                            config.STAT_FISCALYEAR_COL:fiscal_year_renamed},
                            inplace=True)

    # Export to reports
    grant_stats_by_program.to_csv(config.STAT_GRANTS_BY_PROGRAM_FILENAME, 
                                  index=False)
    


def get_shared_projects_by_program_pair(all_grants_data: pd.DataFrame):
    """Get pairs of Key Programs and the count of projects they share"""

    program_field = config.PROGRAM_ID_FIELDNAME
    project_field = config.STAT_CORE_PROJECT_COL

    # Group by core project and combine programs into unordered set
    grouped = all_grants_data.groupby(project_field)[program_field].apply(set)

    # Get combos of programs for each core project with itertools combinations
    all_combinations = grouped.apply(lambda x: list(combinations(
                                        sorted(x), len(x))))

    # Flatten the program combos and count occurrences
    df_shared_programs = (all_combinations.explode()
                          .value_counts().reset_index())

    # Determine the maximum number of programs for a project
    max_programs = (df_shared_programs[program_field]
                        .apply(lambda x: len(x)).max())

    # Create columns dynamically based on the maximum number of programs
    columns = [f'program_{i + 1}' for i in range(max_programs)]
    df_shared_programs[columns] = (pd.DataFrame(
                                    df_shared_programs[program_field].tolist(),
                                    index=df_shared_programs.index))

    # Reorder and rename
    column_names = columns + ['count']
    df_shared_programs = df_shared_programs[column_names]
    df_shared_programs.rename(columns={'count':'unique_core_project_count'},
                              inplace=True)

    # Check if shared programs exist
    if 'program_2' in df_shared_programs:
        # Drop rows where program_2 is NaN or blank to keep only shared projects
        df_shared_programs = df_shared_programs.dropna(subset=['program_2'])
    else:
        # If no shared programs, create a DataFrame with message
        df_shared_programs = pd.DataFrame(
            {'NOTICE': ['No shared projects between programs']})
        
    # Export as report
    shared_programs_filename = config.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME
    df_shared_programs.to_csv(shared_programs_filename, index=False)



def build_program_project_stats(all_grants_data:pd.DataFrame):
    """Create reports with program and project statistics from grants data."""

    print(f"\n---\nPROGRAM/PROJECT STATS:\n"
          f"Generating program and project statistics reports...\n---\n")

    # Define directory to store reports. Create if doesn't already exist
    reports_dir = config.REPORTS_GATHERED_DIR
    if not os.path.exists(reports_dir):
        os.makedirs(reports_dir)  

    # Summary 1: Group by Program and get aggregated grant stats
    get_grant_stats_by_program(all_grants_data)

    # Summary 2: Get pairs of Key Programs that share projects
    get_shared_projects_by_program_pair(all_grants_data)

    print(f"Done! Results can be found in {config.REPORTS_GATHERED_DIR}.")
        



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Load grant.csv as a dataframe
    grant_filepath = config.GRANTS_INTERMED_PATH
    print(f"Generating report statistics on {grant_filepath}...")
    all_cleaned_grants = pd.read_csv(grant_filepath)

    # Run stats module
    build_program_project_stats(all_cleaned_grants)

