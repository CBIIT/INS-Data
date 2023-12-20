"""
main.py
2023-07-26 ZD

Main function for the INS-Data project meant to be run as a single command. 

Inputs required:
- Curated Qualtrics CSV of Key Programs
- iCite database download for publication details
- Reviewed and corrected list of invalid NOFOs/Awards (Optional on subsequent runs)
- All other inputs are gathered via API

Ouputs generated:
- data/01_intermediate/
    - CSVs for all gathered data types, including programs, publications, and 
    projects. These will be passed through data packaging before INS ingestion.
    - Additional CSV outputs and checkpoint files gathered during the process. 
    These are useful for troubleshooting and validation.

- data/02_output/
    - Finalized TSVs for all gathered data types. These are ready for ingestion
    into INS

- reports/:
    - Summary statistic report csvs with high-levelgrant data
    - Report of publications removed from publication.csv with reasons
"""

import os
import sys

import pandas as pd

import config
from modules.data_preparation import load_and_clean_programs
from modules.gather_grants_data import gather_grants_data
from modules.summary_statistics import get_summary_statistics
import modules.gather_publication_data as gpub


def main():
    """Main function for the INS Data Gathering pipeline."""

    # STEP 1: PROGRAMS
    print(f"Loading and processing {config.QUALTRICS_CSV_PATH}...")

    # Load and clean Key Programs
    continue_bool, key_programs_df = load_and_clean_programs(
                                    csv_filepath = config.QUALTRICS_CSV_PATH, 
                                    col_dict = config.QUALTRICS_COLS)

    # If issues are found, user is prompted to stop or continue
    if continue_bool == True:
        # Export cleaned Key Programs file
        key_programs_df.to_csv(config.CLEANED_KEY_PROGRAMS_CSV, index=False)
        print(f"Success! Saved {config.CLEANED_KEY_PROGRAMS_CSV}.")

    else:
        sys.exit("Process ended by user. Cleaned Key Programs csv not saved. "
            "Make manual edits to input CSV and retry.")


    # STEP 2: GRANTS
    
    # Gather, format, and save grants data for each Key Program
    all_cleaned_grants = gather_grants_data(key_programs_df, print_meta=False)


    # STEP 3: STATS
    
    # Run summary statistic module that processes and exports reports as csv
    get_summary_statistics(all_cleaned_grants)
    print(f"Summary statistic reports for grants data successfully generated. \n"
          f"Results can be found in {config.REPORTS_GATHERED_DIR}.\n---")


    # STEP 4: PUBLICATIONS

    # Use existing PMID list if present
    if os.path.exists(config.PROJECT_PMIDS):
        print(f"---\nReusing PMIDs already gathered in {config.PROJECT_PMIDS}.")
        df_pmids = pd.read_csv(config.PROJECT_PMIDS)
    
    else:
        # Gather PMIDs from projects using NIH RePORTER API
        print(f"---\nGathering associated PMIDs for each project from "
              f"NIH RePORTER API...")
        df_pmids = gpub.get_pmids_from_projects(all_cleaned_grants)

        # Store a checkpoint file of all gathered PMIDs and projects
        df_pmids.to_csv(config.PROJECT_PMIDS, index=False)
        print(f"Gathered PMIDs saved to {config.PROJECT_PMIDS}.")

    # Gather publication data and store in chunks of defined size 
    print(f"---\nGathering PubMed data for all PMIDs...\n"
          f"---NOTE: THIS STEP CAN TAKE 8+ HOURS.---\n"
          f"Saving partial files to {config.TEMP_PUBLICATION_DIR} throughout.")
    gpub.build_pmid_info_data_chunks(df_pmids, 
                                chunk_size=config.PUB_DATA_CHUNK_SIZE, 
                                output_folder=config.TEMP_PUBLICATION_DIR)
    
    # Run again to catch any timeouts or problematic PMIDs
    print(f"First pass successful!\n"
          f"Double-checking for any missing publication info...")
    gpub.build_pmid_info_data_chunks(df_pmids, 
                            chunk_size=config.PUB_DATA_CHUNK_SIZE, 
                            output_folder=config.TEMP_PUBLICATION_DIR)
    print(f"Success! PubMed Publication data gathered for all PMIDs.")

    # Gather iCite data for all PMIDs from projects
    print(f"---\nGathering iCite data for all PMIDs...")

    # Use existing iCite-enriched PMIDs if present
    if os.path.exists(config.ICITE_PMID_DATA):
        print(f"---\nReusing iCite data already gathered in {config.ICITE_PMID_DATA}.")
        df_icite = pd.read_csv(config.ICITE_PMID_DATA)
    else:
        df_icite = gpub.get_icite_data_for_pmids(df_pmids,
                                config.ICITE_FILENAME, 
                                config.ICITE_COLUMNS_TO_PULL,
                                chunk_size=250000, chunk_count_est=146)
        # Export iCite data as checkpoint
        df_icite.to_csv(config.ICITE_PMID_DATA, index=False)
        print(f"---\niCite data for PMIDS saved to {config.ICITE_PMID_DATA}.")

    # Load all partial PubMed files back into a single df
    print(f"---\nCombining and merging all publication data...")
    df_pubmed = gpub.load_all_directory_files_to_df(config.TEMP_PUBLICATION_DIR)
    # Reset string-ed dates back to datetime
    df_pubmed['publication_date'] = pd.to_datetime(df_pubmed['publication_date'])

    # Fill in missing PubMed data and columns with iCite
    df_pub_info = gpub.merge_pubmed_icite_pmid_data(df_pubmed, df_icite)
    # Export for troubleshooting
    df_pub_info.to_csv(config.MERGED_PMID_DATA, index=False)
    print(f"Merged PMID data saved to {config.MERGED_PMID_DATA}.")

    # Build final publication df
    print(f"---\nMerging PMID data back to Core Projects...")
    df_publications, df_removed_publications = gpub.merge_and_clean_project_pmid_info(
                                                df_pmids, 
                                                df_pub_info)
    
    # Export final publication data
    df_publications.to_csv(config.PUBLICATIONS_INTERMED_PATH, index=False)
    df_removed_publications.to_csv(config.REMOVED_PUBLICATIONS, index=False)

    print(f"Success! Publication data saved to {config.PUBLICATIONS_INTERMED_PATH}.\n"
          f"Removed publications saved to {config.REMOVED_PUBLICATIONS}")
    print(f"---\n")
    print(f"Total unique Publications saved:  {df_publications['pmid'].nunique():>8}")
    print(f"Total removed Publications:       {len(df_removed_publications):>8}")
    print(f"Project-Publication associations: {len(df_publications):>8}")



if __name__ == "__main__":
    main()
