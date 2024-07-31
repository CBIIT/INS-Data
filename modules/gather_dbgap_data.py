"""
gather_dbgap_data.py
2024-07-26 ZD

This script defines primary function gather_dbgap_data that will accept 
a csv download of studies from the dbGaP advanced search and enrich it with 
data gathered from the dbGaP Study Metadata API and dbGaP SSTR API. The output
will be provided for manual curation of a final dbGaP data and eventual 
integration with a datasets.tsv for INS ingestion.

The output `dbgap_df` and exported dbgap.csv contain columns:
    - TO-DO: list fields here
"""

import os
import sys
import json

import pandas as pd
import requests
from tqdm import tqdm   # for progress bars

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config



def get_dbgap_study_metadata(phs):
    """Fetches study metadata from the dbGaP Study Metadata API in JSON format.

    Args:
        phs: dbGaP phs accession number, e.g. 'phs002790' or 'phs002790.v7.p1'

    Returns:
        A dictionary containing the entire JSON response from the dbGaP API,
            or None if the request fails.
    """ 

    # Base URL for the dbGaP study metadata API
    base_url = "https://submit.ncbi.nlm.nih.gov/dbgap/api/v1/study_config/"

    # Construct the complete URL with the phs accession number
    url = base_url + phs

    # Send a GET request to the API using requests library
    response = requests.get(url)

    # Check for successful response (status code 200)
    if response.status_code == 200:

        # Parse the JSON response using json.loads
        return json.loads(response.text)
    
    else:
        print(f"Error: {phs} - API request failed with status code: "
              f"{response.status_code}")
        return None



def build_raw_study_metadata_records(phs_list:list, export_filename:str):
    """Build a JSON file with raw dbGaP Study Metadata API results for a list  
    of provided phs accessions.  

    Args:
        phs_list (list): List of dbGaP phs accession strings
        export_filename (str): Filename of output json

    Returns:
        None. JSON file is exported.
    """

    print(f"Gathering dbGaP Study Metadata for {len(phs_list)} studies...")

    # Check to see if output file already exists
    if os.path.exists(export_filename):
        print(f"Reusing dbGaP Study Metadata Results "
              f"found in {export_filename}.")

    # Only proceed if file does not already exist
    else:
        
        # Create directory if it does not already exist
        os.makedirs(os.path.dirname(config.DBGAP_META_INTERMED_PATH), 
                    exist_ok=True)

        # Build empty list to fill with results
        all_records = []

        # Iterate through each phs with tqdm progress bar
        for phs in tqdm(phs_list, ncols=80):

            # If already gathered, skip this phs
            if phs in all_records:
                print(f"Duplicate skipped: {phs}")
                continue

            # Pull raw json data from the Study Metadata API
            response = get_dbgap_study_metadata(phs)

            # Proceed if response exists
            if response:
                record = response['data']

                # Add the full phs string used to retrieve the response
                record['full_phs'] = phs

                # Add record to the combined list
                all_records.append(record)

            # TO DO: Gather any failed phs into a report for reference

        # Export as json
        with open(export_filename, 'w') as outfile:
            json.dump(all_records, outfile, indent=2)

        print(f"---\nSuccess! Study Metadata for {len(all_records)} "
              f"saved to {export_filename}.")



def get_principle_investigators(record):
    """Parse a dbGaP Study Metadata record to get all Principal Investigators 
    as a semicolon-separated list-like string.

    Args:
        record (str): Study Metadata API response for a single study.
            Must contain 'attribute' field.

    Returns:
        str: Semicolon-separated list of Principal Investigators
    """

    # Build empty list to fill with PIs
    pi_list = []

    # List of all possible 'title' strings of interest
    title_list = ['Principal Investigator', 
                  'Principal investigator',
                  'Principal Investigators',
                  'Principal investigators',]

    # Iterate through all attribution records to find matching titles
    for attribution in record['attribution']:
        if attribution.get('title') in title_list:

            # Collect the name field and add to running list
            pi_list.append(attribution['name'])

    # Connect all results as strings with separator
    pi_output = ';'.join(pi_list)

    return pi_output



def get_cited_publications(record):
    """Parse a dbGaP Study Metadata record to get all cited publications 
    as a semicolon-separated list-like string. 

    Args:
        record (str): Study Metadata API response for a single study. 
            Must contain 'reference' field.

    Returns:
        str: Semicolon-separated list of publications as pmids or titles
    """

    # Build empty list to fill with citations
    citation_list = []

    # Iterate through all references listed
    for ref in record['reference']:

        # Gather any PMIDs available 
        if 'pmid' in ref:
            citation_list.append(ref['pmid'])
        
        # If PMID not available, collect the article title
        elif 'title' in ref:
            citation_list.append(ref['title'])

    # Connect all results as strings with separator
    citation_output = ';'.join(map(str, citation_list))

    return citation_output



def get_funding_attributions(record):
    """Parse a dbGaP Study Metadata record to get all funding information 
    as a semicolon-separated list-like string. 

    Args:
        record (str): Study Metadata API response for a single study. 
            Must contain 'attribution' field.

    Returns:
        str: Semicolon-separated list of provided funding information
    """

    # Build empty list to fill with funding
    funding_list = []

    # List of all possible 'title' strings of interest
    title_list = ['Funding Source']

    # Iterate through all attribution records to find matching titles
    for attribution in record['attribution']:
        if attribution.get('title') in title_list:

            # Collect the name field and add to running list
            funding_list.append(attribution['name'])

    # Connect all results as strings with separator
    funding_output = ';'.join(map(str, funding_list))

    return funding_output



def get_external_study_urls(record):
    """Parse a dbGaP Study Metadata record to get all provided study urls
    as a semicolon-separated list-like string. 

    Args:
        record (str): Study Metadata API response for a single study. 
            Must contain 'study_url' field.

    Returns:
        str: Semicolon-separated list of provided study urls
    """

    # Build empty list to fill with urls
    url_list = []

    # Iterate through all study_url fields
    for link in record['study_url']:

        # Collect the url field and add to running list
        url = link.get('url', None)
        url_list.append(url)

    # Connect all results as strings with separator
    url_output = ';'.join(map(str, url_list))

    return url_output



def join_list(items, sep=';'):
    """Joins a list of items into a delimited string.

    Args:
        items (list): List of items to join
        sep (str): Character to separate items. Default semicolon

    Returns:
        str: Comma-separated string of items
    """

    # Build empty list to fill with items
    item_list = []

    # Iterate through items and add to running list
    for item in items:
        item_list.append(item)
    
    # Map values to string and connect with separator
    string_list = sep.join(map(str, item_list))

    return string_list



def clean_dbgap_study_metadata(record:str) -> dict:
    """Process study metadata json into a standardized, flat dictionary. 

    Args:
        record (str): Nested JSON results from the dbGaP Study Metadata API

    Returns:
        study_metadata_info (dict): Flat dict of selected metadata
    """

    study_metadata = {
        'full_phs': record['full_phs'],
        'principal_investigator': get_principle_investigators(record),
        'cited_publications': get_cited_publications(record),
        'funding_source': get_funding_attributions(record),
        'study_type': join_list(record['study_type']),
        'external_study_url': get_external_study_urls(record),
        'gene_keywords': join_list(record['gene']),
        'disease_keywords': join_list(record['disease']),
    }

    return study_metadata



def get_dbgap_sstr_metadata(phs:str, api_option:str='summary'):
    """Fetches study metadata from the dbGaP SSTR API in JSON format.

    Args:
        phs (str): dbGaP phs accession number, e.g. 'phs002790.v7.p1' 
            If base accession is used e.g. ('phs002790'), the response will use
            the most recent available version. 
        api_option (str): 'subjects' or 'summary'. Determines which SSTR API 
            response type will be called. Default 'summary'.

    Returns:
        A dictionary containing the entire JSON response from the dbGaP API,
            or None if the request fails.
    """ 

    # Base URL for the dbGaP study metadata API
    base_url = "https://www.ncbi.nlm.nih.gov/gap/sstr/api/v1/study/"

    # Construct the complete URL with the phs accession number and API type
    url = base_url + phs + "/" + api_option

    # Send a GET request to the API using requests library
    response = requests.get(url)

    # Check for successful response (status code 200)
    if response.status_code == 200:

        # Parse the JSON response using json.loads
        return json.loads(response.text)
    
    else:
        print(f"Error: {phs} - API request failed with status code: "
              f"{response.status_code}")
        return None



def build_raw_sstr_records(phs_list:list, export_filename:str):
    """Build a JSON file with raw dbGaP SSTR API results for a list  
    of provided phs accessions.  

    Args:
        phs_list (list): List of dbGaP phs accession strings
        export_filename (str): Filename of output json

    Returns:
        None. JSON file is exported.
    """

    print(f"Gathering dbGaP SSTR Metadata for {len(phs_list)} studies...")

    # Check to see if output file already exists
    if os.path.exists(export_filename):
        print(f"Reusing dbGaP SSTR Results "
              f"found in {export_filename}.")

    # Only proceed if file does not already exist
    else:
        
        # Create directory if it does not already exist
        os.makedirs(os.path.dirname(config.DBGAP_SSTR_INTERMED_PATH), 
                    exist_ok=True)

        # Build empty list to fill with results
        all_records = []

        # Iterate through each phs # ADD TQDM PROGRESS BAR
        for phs in phs_list:

            # If already gathered, skip this phs
            if phs in all_records:
                print(f"Duplicate skipped: {phs}")
                continue

            # Pull raw json data from the SSTR API
            record = get_dbgap_sstr_metadata(phs)

            # Proceed if response exists
            if record:

                # Add the full phs string used to retrieve the response
                record['full_phs'] = phs

                # Add record to the combined list
                all_records.append(record)

            # TO DO: Gather any failed phs into a report for reference

        # Export as json
        with open(export_filename, 'w') as outfile:
            json.dump(all_records, outfile, indent=2)

        print(f"---\nSuccess! SSTR Metdata for {len(all_records)} "
              f"saved to {export_filename}.")



def get_sstr_count(record:str, count_type:str):
    """Sum count of participants (aka subjects) or samples across all SSTR 
    study consent groups for a single dbGaP study.

    Args:
        record (str): JSON SSTR API response for a single study
        count_type (str): 'participant' or 'sample'. Value to sum and return
    """

    # Derive targeted fields from count_type
    if count_type == 'participant':
        count_field = 'subject_count'
    elif count_type == 'sample':
        count_field = 'sample_count'

    # Throw an error if incorrect count_type is used
    else:
        raise ValueError(f"Invalid 'count_type': '{count_type}'. Must be "
                         f"either 'participant' or 'sample'.")
        

    # Find all consent groupings listed in the response
    consent_groups = record['study']['consent_groups']

    # Start a sum counter for across groups
    total_count = 0
    
    # Pull and sum the defined type
    for group in consent_groups:
        total_count += group[count_field]
        
    return total_count



def get_consent_codes(record:str):
    """Pull all consent code values present within a single dbGaP study.

    Args:
        record (str): JSON SSTR API response for a single study
    """

    # Build empty set to fill with consent codes
    code_list = set()

    # Find all consent groupings listed in the response
    consent_groups = record['study']['consent_groups']

    # Find all consent codes and add code short name to running list
    for group in consent_groups:
        code_list.add(group['short_name'])

    # Connect all results as strings with separator
    consent_output = ';'.join(code_list)

    return consent_output



def get_assay_method(record:str):
    """Pull all sample assay values present within a single dbGaP study.

    Args:
        record (str): JSON SSTR API response for a single study
    """

    # Build empty set to fill with assay methods
    method_list = set()

    # Find all consent and sample groupings listed in the response
    sample_groups = record['study_stats']['cnt_by_consent_and_sample_use']

    # Find all assay methods and add to running list
    for group in sample_groups:
        method_list.add(group['sample_use'])

    # Connect all results as strings with separator
    method_output = ';'.join(method_list)

    return method_output



def clean_dbgap_sstr_metadata(record:str) -> dict:
    """Process SSTR response JSON into a standardized, flat dictionary. 

    Args:
        record (str): JSON results from the dbGaP SSTR API

    Returns:
        sstr_metadata: (dict): Flat dict of selected metadata
    """

    sstr_metadata = {
        'full_phs': record['full_phs'],
        'participant_count': get_sstr_count(record, count_type='participant'),
        'sample_count': get_sstr_count(record, count_type='sample'),
        'limitations_for_reuse': get_consent_codes(record),
        'assay_method': get_assay_method(record)
    }

    return sstr_metadata



def gather_dbgap_data(input_csv:str):
    """Gather dbGaP study data from multiple sources.

    Args:
        input_csv (str): Filepath of CSV containing download of studies of 
            interest from the dbGaP Advanced Search UI

    Returns:
        None. Exports a csv of dbGaP datasets
    """

    print(f"\n---\nDATASETS - dbGaP:\n"
        f"Gathering, formatting, and saving dbGaP dataset metadata...\n---\n")

    # Read the dbGaP CSV download
    df = pd.read_csv(input_csv)

    # Define list of fields from the CSV download to keep in results
    csv_cols_to_keep = [
        'accession',
        'name',
        'description',
        'Study Disease/Focus',
        #'Study Molecular Data Type',  # same values as SSTR assay method
        'Release Date',
        'Ancestry (computed)',
        'Related Terms'
    ]

    # Get list of phs accessions of interest 
    phs_list = df['accession']

    #---
    # STUDY METADATA API
    print(f"\nProcessing dbGaP Study Metadata API results...")

    # Get raw responses from the dbGaP Study Metadata API and store as file
    build_raw_study_metadata_records(phs_list, config.DBGAP_META_INTERMED_PATH)

    # Load the Study Metadata records from stored file
    with open(config.DBGAP_META_INTERMED_PATH, 'r') as file:
        study_metadata_records = json.load(file)

    # Build emtpy list to fill with processed records
    processed_meta_record_list = []

    # Iterate through study records for processing and add to running list
    for meta_record in tqdm(study_metadata_records, ncols=80):
        processed_record = clean_dbgap_study_metadata(meta_record)
        processed_meta_record_list.append(processed_record)

    # Convert to dataframe
    metadata_df = pd.DataFrame(processed_meta_record_list)

    # Merge dbGaP downloaded CSV results with Study Metadata results
    meta_merged_df = pd.merge(df[csv_cols_to_keep], metadata_df, 
                         left_on='accession', right_on='full_phs', how='left')
    meta_merged_df.drop('full_phs', axis=1, inplace=True)

    #---
    # SSTR API

    # Get raw responses from the dbGaP SSTR API and store as file
    build_raw_sstr_records(phs_list, config.DBGAP_SSTR_INTERMED_PATH)

    # Load the SSTR records from stored file
    with open(config.DBGAP_SSTR_INTERMED_PATH, 'r') as file:
        study_sstr_records = json.load(file)

    # Build emtpy list to fill with processed records
    processed_sstr_list = []

    # Iterate through sstr records for processing and add to running list
    for sstr_record in tqdm(study_sstr_records, ncols=80):
        processed_sstr = clean_dbgap_sstr_metadata(sstr_record)
        processed_sstr_list.append(processed_sstr)

    # Convert to dataframe
    sstr_df = pd.DataFrame(processed_sstr_list)

    # Merge SSTR fields with CSV download + Study Metadata
    merged_df = pd.merge(meta_merged_df, sstr_df,
                         left_on='accession', right_on='full_phs', how='left')
    merged_df.drop('full_phs', axis=1, inplace=True)

    # Export as CSV 
    merged_df.to_csv(config.DBGAP_PROCESSED_PATH, index=False)

    print(f"\nSuccess! dbGaP data saved to {config.DBGAP_PROCESSED_PATH}.\n")



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Gather, process, and save dbgap dataset metadata
    gather_dbgap_data(config.DBGAP_INPUT_CSV)
