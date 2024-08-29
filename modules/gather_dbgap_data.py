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
import re

import pandas as pd
import requests
from tqdm import tqdm   # for progress bars
from Bio import Entrez  # for PubMed API

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config



def get_dbgap_api_data(phs:str, api_type:str):
    """Fetches data from the dbGaP Study Metadata or SSTR API in JSON format.

    Args:
        phs (str): dbGaP phs accession number, e.g. 'phs002790' or 'phs002790.v7.p1'
        api_type (str): 'study_metadata', 'sstr_summary' or 'sstr_subjects'. 
            Determines which dbGaP API will be called.

    Returns:
        A dictionary containing the entire JSON response from the dbGaP API,
            or None if the request fails.
    """ 

    if api_type == 'study_metadata':
        # Base URL for the dbGaP study metadata API
        base_url = "https://submit.ncbi.nlm.nih.gov/dbgap/api/v1/study_config/"
        url = base_url + phs

    elif api_type == 'sstr_summary':
        # Base URL for the dbGaP SSTR Summary API
        base_url = "https://www.ncbi.nlm.nih.gov/gap/sstr/api/v1/study/"
        url = base_url + phs + "/summary"

    elif api_type == 'sstr_subjects':
        # Base URL for the dbGaP SSTR Subjects (aka participants) API
        base_url = "https://www.ncbi.nlm.nih.gov/gap/sstr/api/v1/study/"
        url = base_url + phs + "/subjects"

    else:
        raise ValueError(f"Invalid 'api_type': '{api_type}'. Must be "
                         f"either 'study_metadata', 'sstr_summary', "
                         f"or 'sstr_subjects'.")

    # Call the API with requests
    response = requests.get(url)

    # Parse the response as JSON
    record = json.loads(response.text)

    # Check for unsuccessful response (status code 200 is success)
    if response.status_code != 200:
        
        tqdm.write(f"Error: {phs} - API request failed with status code: "
              f"{response.status_code}")
        
        # Add a field listing the response code for errors
        record['error']['response_code'] = response.status_code

    return record



def build_bulk_raw_dbgap_api_data(phs_list:list, 
                                  api_type: str,
                                  export_filename:str, 
                                  error_report_filename:str):
    """Build a JSON file with unprocessed dbGaP API results for a list  
    of provided phs accessions.  

    Args:
        phs_list (list): List of dbGaP phs accession strings
        api_type (str): 'study_metadata', 'sstr_summary' or 'sstr_subjects'. 
            Determines which dbGaP API will be called.
        export_filename (str): Filename of output json
        error_report_filename (str): filename of output csv with api errors

    Returns:
        None. JSON file is exported.
    """

    # Check to see if output file already exists
    if os.path.exists(export_filename):
        print(f"Reusing dbGaP API bulk results for {api_type}"
              f"found in {export_filename}.")

    # Only proceed if file does not already exist
    else:
        
        # Create directory if it does not already exist
        os.makedirs(os.path.dirname(export_filename), 
                    exist_ok=True)

        # Build empty lists to fill with results
        success_records = []
        error_records = []

        # Iterate through each phs with tqdm progress bar
        for phs in tqdm(phs_list, ncols=80):

            # If already gathered, skip this phs
            if phs in success_records or phs in error_records:
                tqdm.write(f"Duplicate skipped: {phs}")
                continue

            # Pull raw json data from the API
            response = get_dbgap_api_data(phs, api_type)

            # Proceed only if response exists
            if response:

                # Handle successful responses
                # Study Metadata API uses 'data' as top-level field
                if response.get('data'):
                    record = response['data']
                    record_type = 'success'

                # SSTR API uses 'study' and 'study_stats' as top-level fields
                elif response.get('study') or response.get('study_stats'):
                    record = response
                    record_type = 'success'
                
                # Handle error response
                elif response.get('error'):
                    record = response['error']
                    record_type = 'error'

                else:
                    raise KeyError(f"Expected top-level field not found in API "
                                   f"response for {phs}.")

                # Add the full phs string used to retrieve the response
                record['full_phs'] = phs

                # Add record to appropriate running list
                if record_type == 'success':
                    success_records.append(record)
                if record_type == 'error':
                    error_records.append(record)

                # Throw warning if there is no response at all
            else: 
                raise Warning(f"Error: {phs} - No API response.")

        # Export error report as csv via pandas
        if len(error_records)>0:

            # Create directory if it does not already exist
            os.makedirs(os.path.dirname(error_report_filename), 
                    exist_ok=True)
            
            # Load error list to pandas dataframe
            error_df = pd.DataFrame(error_records)
            
            # Move phs to first column
            error_df.insert(0, 'full_phs', error_df.pop('full_phs'))

            # Export to csv
            error_df.to_csv(error_report_filename, index=False)

            print(f"Details for {len(error_records)} studies with API failures "
                  f"saved to {error_report_filename}.\n")

        # Only export results if any exist
        if len(success_records) > 0:

            # Export good data as json
            with open(export_filename, 'w') as outfile:
                json.dump(success_records, outfile, indent=2)

            print(f"---\nSuccess! dbGaP {api_type} data for "
                  f"{len(success_records)} studies saved "
                  f"to {export_filename}.\n")
        else: 
            raise Warning(f"No studies found using API: {api_type}.")



def detect_attribute_fieldnames(raw_records, title_type):
    """Parse dbGaP Study Metadata records to get all attribute fieldnames likely
    to include information about a specified field of interest.
    Principal Investigators (PIs) and Funding Sources are supported.

    Args:
        raw_records (str): JSON Study Metadata API responses for all studies.
            Records must contain 'attribute' field.
        title_type (str): 'pi' or 'funding'. Specifies which sets of regex 
            patterns to use when detecting relevant titles. 

    Returns:
        set: Set of all relevant titles within record[`attribute`]
    """

    # Define regex patterns to detect relevant attribution fieldnames
    # Structure as dict with re.compile flags as values

    if title_type == 'pi':
    # Define Principal Investigator search patterns
        patterns = {
            r"(?s:.)*PI(?s:.)*": re.NOFLAG,
            r"(?s:.)*Prin(?s:.)*In(?s:.)*": re.IGNORECASE,
            r"(?s:.)*Lead(?s:.)*In(?s:.)*": re.IGNORECASE
        }

        # Define regex patterns to exclude from detection
        exclusion_patterns = {
            r".*for princip*": re.IGNORECASE,
            r".*of princip*": re.IGNORECASE,
        }
    
    elif title_type == 'funding':
    # Define Funding search patterns
        patterns = {
            r"(?s:.)*Funding(?s:.)*": re.IGNORECASE,
            r"(?s:.)*Grant(?s:.)*": re.IGNORECASE,
            r"(?s:.)*Award(?s:.)*": re.IGNORECASE,
        }

        # Define regex patterns to exclude from detection
        exclusion_patterns = {}

    else:
        raise ValueError(f"Invalid 'title_type': '{title_type}'. Must be "
                         f"either 'pi' or 'funding'.")

    # Build empty set to fill with PI-like fieldname titles
    title_set = set()

    # Iterate through each study record in the input json records
    for record in raw_records:

        # Iterate through all attribution titles within single study record
        for attribution in record['attribution']:
            title = attribution.get('title')

            # Continue if any titles exist
            if title:

                # Iterate through the detection patterns
                for pattern, flag in patterns.items():
                    
                    # Check if the given title matches the given pattern
                    compiled_pattern = re.compile(pattern, flags=flag)

                    # Check for exclusions after a pattern is matched
                    if compiled_pattern.search(title): 
                        is_excluded = any(re.search(exclusion, title, flag) 
                                          for exclusion, flag 
                                          in exclusion_patterns.items())

                        if not is_excluded:
                            # Add to running list to title fieldnames
                            title_set.add(title)

    return title_set



def get_principal_investigators(record, pi_title_set):
    """Parse a dbGaP Study Metadata record to get all Principal Investigators 
    as a semicolon-separated list-like string.

    Args:
        record (str): Study Metadata API response for a single study.
            Must contain 'attribute' field.
        pi_title_set (set): Set of attribution titles likely to include names
            of Principal Investigators (PIs).

    Returns:
        str: Semicolon-separated list of Principal Investigators
    """

    # Build empty list to fill with PIs
    pi_list = []

    # Customization option to add or exclude titles to the regex-generated set
    add_titles = set(['Principal Investigator',
                     'Senior Principal Scientist',
                     'John Prensner, MD PhD'])
    remove_titles = set(['Known bad title',])

    pi_title_set.update(add_titles)
    pi_title_set = pi_title_set - remove_titles

    # Convert set to list for downstream iteration
    pi_title_list = list(pi_title_set)

    # Iterate through all attribution records to find matching titles
    for attribution in record['attribution']:
        if attribution.get('title') in pi_title_list:

            # Collect the name field and add to running list
            pi_list.append(attribution['name'])

    # Connect all results as strings with separator
    pi_output = ';'.join(pi_list)

    return pi_output



def get_pmid_from_reference(search_string, type):
    """Retrieves the PubMed ID (PMID) for a given title or author list.

    Args:
        search_term (str): The title of the article
        type (str): 'title' or 'authors'. Input search type

    Returns:
        int: PMID of first resulting article, or None if no results are found.

    Raises:
        ValueError: If the input is not a string.
        RuntimeError: If there's an unexpected error during the search.
    """

    # Get user email from hidden local env file. Use default if not defined
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')

    # Check for input type validity
    if not isinstance(search_string, str):
        raise ValueError("Article title must be a string.")

    # Remove special character formatted words to avoid exact search issues
    words = search_string.split()
    filtered_words = []
    for word in words:
        if "&#" not in word:
            filtered_words.append(word)
    filtered_string = " ".join(filtered_words)
    
    try:
        # Call the PubMed e-utilities via Entrez
        handle = Entrez.esearch(db="pubmed", term=filtered_string, 
                                retmax=1, field=type)
        record = Entrez.read(handle)
        handle.close()

        # Return none if no matching records found
        if record["Count"] == "0":
            return None

        # Get the first resulting PMID as an integer
        pmid_out = int(record["IdList"][0])

        return pmid_out

    # Throw an error and identify any Entrez issues
    except Exception as e:
        tqdm.write(f"Error retriving PMID for {type} {filtered_string}: {str(e)}")
        return None



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

        # Set default pmid to None to avoid value error
        pmid = None

        # Gather any PMIDs directly available 
        if 'pmid' in ref:
            citation_list.append(ref['pmid'])
        
        # If not directly availabe, use Entrez e-utilities to get PMID
        else:
            # First try to find a PMID from the article title
            if 'title' in ref:
                pmid = get_pmid_from_reference(ref['title'], 'title')

                # Then try to find a PMID from the authors
                if pmid == None and 'authors' in ref:
                    pmid = get_pmid_from_reference(ref['authors'], 'authors')

            # Use PMID if found or article title if not
            if pmid: 
                citation_list.append(pmid)
            else: 
                citation_list.append(f"No PMID: {ref['title']}")

    # Connect all results as strings with separator
    citation_output = ';'.join(map(str, citation_list))

    return citation_output



def get_funding_attributions(record, funding_title_set):
    """Parse a dbGaP Study Metadata record to get all funding information 
    as a semicolon-separated list-like string. 

    Args:
        record (str): Study Metadata API response for a single study. 
            Must contain 'attribution' field.
        funding_title_set (set): Set of attribution titles likely to include
            details about study funding.

    Returns:
        str: Semicolon-separated list of provided funding information
    """

    # Build empty list to fill with funding
    funding_list = []

    # Customization option to add or exclude titles to the regex-generated set
    add_titles = set(['Funding Source',
                      'Founding Source',])
    remove_titles = set(['Known bad title',])

    funding_title_set.update(add_titles)
    funding_title_set = funding_title_set - remove_titles

    # Convert set to list for downstream iteration
    funding_title_list = list(funding_title_set)

    # Iterate through all attribution records to find matching titles
    for attribution in record['attribution']:
        if attribution.get('title') in funding_title_list:

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



def clean_dbgap_study_metadata(record:str, 
                               pi_title_set:set, 
                               funding_title_set:set) -> dict:
    """Process study metadata json into a standardized, flat dictionary. 

    Args:
        record (str): Nested JSON results from the dbGaP Study Metadata API

    Returns:
        study_metadata_info (dict): Flat dict of selected metadata
    """

    study_metadata = {
        'full_phs': record['full_phs'],
        'principal_investigator': get_principal_investigators(record, 
                                                              pi_title_set),
        'cited_publications': get_cited_publications(record),
        'funding_source': get_funding_attributions(record, 
                                                   funding_title_set),
        'study_type': join_list(record['study_type']),
        'external_study_url': get_external_study_urls(record),
        'gene_keywords': join_list(record['gene']),
        'disease_keywords': join_list(record['disease']),
    }

    return study_metadata



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



def get_dbgap_url(accession:str):
    """Build a url to the dbGaP study page using phs accession."""

    base_url = 'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id='

    # Check that phs accession string begins with phs
    if not accession.startswith('phs'):
        raise ValueError(f"Invalid phs accession: '{accession}'. Must begin "
                         f"with 'phs'.")

    else:
        # Get the shortened accession without versioning
        short_accession = accession[0:9]
        
        # Combine with base url for study page
        url = base_url + short_accession

        return url



def get_gpa_and_doc(dbgap_df:pd.DataFrame):
    """Add Grant Program Administrator (GPA) names and NCI Division, Offices, 
    or Centers (DOC) labels for each dbGaP dataset. 

    Args:
        dbgap_df (pd.DataFrame): pandas dataframe containing dbGaP information.
            Must contain 'accession' column with dbGaP phs accessions.
    """

    # Load GPA table received from ODS
    # This hold GPA names and associated phs accessions
    gpa_df = pd.read_csv(config.DBGAP_GPA_LIST)

    # Load gpa/doc lookup table
    # This holds GPA names and associated DOCs
    gpa_doc_lut_df = pd.read_csv(config.DBGAP_GPA_DOC_LUT)

    # Merge the two GPA tables into one lookup table
    merged_lut_df = pd.merge(left=gpa_df[['Accession','Primary GPA']], 
                             right=gpa_doc_lut_df[['gpa', 'doc']],
                             how='left', 
                             left_on='Primary GPA', right_on='gpa')
    
    # Add the GPA and DOC columns to the main df
    df_out = pd.merge(left=dbgap_df, 
                      right=merged_lut_df[['Accession','gpa','doc']],
                      how='left', 
                      left_on='accession', right_on='Accession')
    
    # Drop the extra column
    df_out.drop(columns='Accession',inplace=True)

    return df_out



def select_and_rename_columns(df:pd.DataFrame):
    """Standardize columns for the datasets output file. 
    NOTE: Consider integrating this into `config.py` and the 
    `package_output_data.py` in the future for consistency.

    Args:
        df (pd.DataFrame): pandas dataframe of merged dataset metadata
    """

    # Define columns to keep, reorder, and rename
    # Column names should match data model properties
    dataset_cols = {
        'name': 'dataset_title',
        'description': 'description',
        'accession': 'dbGaP_phs',
        'dbGaP_URL': 'dbGaP_URL',
        'principal_investigator': 'PI_name',
        'gpa': 'GPA',
        'doc': 'dataset_doc',
        'cited_publications': 'dataset_pmid',
        'funding_source': 'funding_source',
        'Release Date': 'release_date',
        'limitations_for_reuse': 'limitations_for_reuse',
        'assay_method': 'assay_method',
        'study_type': 'study_type',
        'Study Disease/Focus': 'primary_disease',
        'participant_count': 'participant_count',    
        'sample_count': 'sample_count',
        'external_study_url': 'study_links',
        'gene_keywords': 'related_genes',
        'disease_keywords': 'related_diseases',
        'Related Terms': 'related_terms',
    }

    # Rename and reorder columns
    df = df.rename(columns=dataset_cols)

    # Keep only defined columns and drop others
    df = df[list(dataset_cols.values())]

    return df



def build_dbgap_df_from_json(api_type: str,
                             input_filepath: str):
    """Process a JSON file of bulk raw dbGaP API results into a dataframe
    with a subset of information relevant for INS. 

    Args:
        api_type (str): 'study_metadata', or 'sstr_summary'. 
            Determines which processing function will be run. 
            'sstr_subjects' is not yet used or supported.
        input_filepath (str): Filepath of JSON file containing API results

    """

    # Load the API response records from stored file
    with open(input_filepath, 'r') as file:
        raw_records = json.load(file)

    # Build empty list to fill with processed records
    processed_record_list = []

    # Detect attribute titles for Study Metadata fields of interest
    if api_type == 'study_metadata':
        pi_title_set = detect_attribute_fieldnames(raw_records, 'pi')
        funding_title_set = detect_attribute_fieldnames(raw_records, 'funding')

    # Iterate through records for processing and add to running list
    for record in raw_records:

        # Use the correct functions to clean and process data
        if api_type == 'sstr_summary':
            processed_record = clean_dbgap_sstr_metadata(record)

        elif api_type == 'study_metadata':
            processed_record = clean_dbgap_study_metadata(record,
                                                          pi_title_set,
                                                          funding_title_set)

        else:
            raise ValueError(f"Invalid 'api_type': '{api_type}'. Must be "
                            f"either 'study_metadata', or 'sstr_summary'. "
                            f" 'sstr_subjects' is not yet supported.")

        processed_record_list.append(processed_record)

    # Convert to dataframe
    df = pd.DataFrame(processed_record_list)

    return df



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
        'Release Date',
        'Ancestry (computed)',
        'Related Terms'
    ]

    # Get list of phs accessions of interest 
    phs_list = df['accession']


    print(f"\nProcessing dbGaP Study Metadata API results...\n")

    # Get raw responses from the dbGaP Study Metadata API and store as file
    build_bulk_raw_dbgap_api_data(phs_list, 
                                  'study_metadata',
                                  config.DBGAP_META_INTERMED_PATH,
                                  config.DBGAP_META_ERRORS)

    # Process Study Metadata
    metadata_df = build_dbgap_df_from_json('study_metadata',
                                           config.DBGAP_META_INTERMED_PATH)


    print(f"\nProcessing dbGaP SSTR API results...\n")

    # Get raw responses from the dbGaP SSTR API and store as file
    build_bulk_raw_dbgap_api_data(phs_list, 
                                  'sstr_summary',
                                  config.DBGAP_SSTR_INTERMED_PATH,
                                  config.DBGAP_SSTR_ERRORS)
    
    sstr_df = build_dbgap_df_from_json('sstr_summary',
                                       config.DBGAP_SSTR_INTERMED_PATH)


    print(f"\nMerging dbGaP data from multiple sources...")

    # Merge CSV download + Study Metadata
    meta_merged_df = pd.merge(df[csv_cols_to_keep], metadata_df, 
                         left_on='accession', right_on='full_phs', how='left')
    meta_merged_df.drop('full_phs', axis=1, inplace=True)

    # Merge SSTR fields with CSV download + Study Metadata
    merged_df = pd.merge(meta_merged_df, sstr_df,
                         left_on='accession', right_on='full_phs', how='left')
    merged_df.drop('full_phs', axis=1, inplace=True)


    print(f"\nFormatting dbGaP dataset output...")

    # Add dbGaP URL column
    merged_df['dbGaP_URL'] = merged_df['accession'].apply(get_dbgap_url)

    # Add GPA and DOC columns
    dbgap_df = get_gpa_and_doc(merged_df)

    # Reorder, rename, and select output columns
    dbgap_df = select_and_rename_columns(dbgap_df)

    # Export final merged df as CSV 
    dbgap_df.to_csv(config.DBGAP_PROCESSED_PATH, index=False)

    print(f"\nSuccess! dbGaP data saved to {config.DBGAP_PROCESSED_PATH}.\n")



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Gather, process, and save dbgap dataset metadata
    gather_dbgap_data(config.DBGAP_INPUT_CSV)
