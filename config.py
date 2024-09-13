"""
config.py for INS Data Processing
2023-07-26 ZD
Make changes here to affect variables throughout repo
"""

from datetime import datetime

# FILEPATH CONFIGURATION

# Edit version and type below for each new Qualtrics file received
# Inputs and outputs will use this versioning
# Version must match suffix in input filename

QUALTRICS_VERSION = "2024-04-24"    # <-- CHANGE VERSION HERE
QUALTRICS_TYPE = "manual_fix"       # <-- Define "raw" or "manual_fix" type of the input csv

# Version of bulk download from iCite
ICITE_VERSION = "2024-03"           # <-- CHANGE VERSION HERE

# Version of dbGaP seearch results download (download date)
DBGAP_CSV_VERSION = "2024-08-30"   # <-- CHANGE VERSION HERE

# An override date can be used instead of today's date for pulling and saving data versions
# This is useful when running downstream modules on grants data gathered before today

OVERRIDE_DATE = None    # <-- Optional. Define override date (e.g. "2023-12-14"). Default None.



# --- DO NOT EDIT BELOW FOR ROUTINE DATA GATHERING ---

INPUT_DIR = "data/00_input/"
INTERMED_DIR = "data/01_intermediate/"
OUTPUT_DIR = "data/02_output/"



# Qualtrics Programs input path
QUALTRICS_CSV_PATH = INPUT_DIR + "qualtrics/qualtrics_output_" + QUALTRICS_VERSION +"_"+ QUALTRICS_TYPE + ".csv"

# Add timestamp to note when grants were gathered from API
# The same Qualtrics input file can have different outputs depending upon API gathering date
TIMESTAMP = "gathered-"+datetime.now().strftime("%Y-%m-%d")

# Use optional override date if provided
if OVERRIDE_DATE:
    TIMESTAMP = "gathered-" + OVERRIDE_DATE
    print(f"\n---TIMESTAMP OVERRIDE IN USE---\n"
          f"---Performing action using {OVERRIDE_DATE} instead of current timestamp.---\n"
          f"---Change the OVERRIDE_DATE in config.py to None to restore default behavior.---\n\n")

# Versioned directories for intermediates and outputs
GATHERED_DIR = INTERMED_DIR + QUALTRICS_VERSION +"/"+ TIMESTAMP
OUTPUT_QUALTRICS_DIR = OUTPUT_DIR + QUALTRICS_VERSION
OUTPUT_GATHERED_DIR = OUTPUT_QUALTRICS_DIR +"/"+ TIMESTAMP

# Reports directories
REPORTS_DIR = "reports/" + QUALTRICS_VERSION
REPORTS_GATHERED_DIR = REPORTS_DIR +"/"+ TIMESTAMP

# Data validation for QA. Filename is tagged with the date it was generated
DATA_VALIDATION_EXCEL = REPORTS_GATHERED_DIR +"/"+ "INS_DataValidation_Generated_" + datetime.now().strftime(("%Y-%m-%d"))+".xlsx"

# Programs output
PROGRAMS_INTERMED_PATH = INTERMED_DIR + QUALTRICS_VERSION +"/"+ "program.csv"
PROGRAMS_OUTPUT_PATH = OUTPUT_GATHERED_DIR +"/"+ "program.tsv"

# Grants output
GRANTS_INTERMED_PATH = GATHERED_DIR +"/"+ "grant.csv"
GRANTS_OUTPUT_PATH = OUTPUT_GATHERED_DIR +"/"+ "grant.tsv"

# Projects output
PROJECTS_INTERMED_PATH = GATHERED_DIR +"/"+ "project.csv"
PROJECTS_OUTPUT_PATH = OUTPUT_GATHERED_DIR +"/"+ "project.tsv"



# ---
# DATA PREPARATION CONFIGURATION

# Dictionary of old:new column names to keep from Qualtrics
QUALTRICS_COLS = {
    "Name of Key Program": "program_name",
    "Acronym for key program": "program_acronym",
    "Focus Area (select all that apply)": "focus_area",
    "DOC": "doc",
    "Primary Contact (PI)": "contact_pi",
    "Primary Contact (PI) email": "contact_pi_email",
    "NIH Contact (Program Officer/Program Director)": "contact_nih",
    "NIH Contact (Program Officer/Program Director) email": "contact_nih_email",
    "NOFO number (eg. format as \"RFA-CA-00-000\") (If more than one, separate with ; semicolon)": "nofo",
    "Grant/Award number {parent award FORMAT LL#CA######, eg. UG3CA260607} (If more than one, separate with ; semicolon)": "award",
    "Link to program website": "program_link",
    "Link to data or DCC if available": "data_link",
    "What type of cancer is the primary focus of the program? (Check all that\napply)": "cancer_type",
    "Login ID": "login_id"
}

# Generic value to use when no specific cancer type is specified
PROGRAM_FILLER_CANCER_TYPE = 'Broad Cancer Types'

# Dictionary of specific old:new values to replace within data
PROGRAM_VALUE_REPLACEMENTS = {"This program focuses on cancer broadly - not limited to a primary cancer type": PROGRAM_FILLER_CANCER_TYPE}

# Dictionary of column_name:filler_value to replace blank values within specific columns
PROGRAM_BLANK_REPLACEMENTS ={
    'focus_area': 'No Focus Area',
    'cancer_type': PROGRAM_FILLER_CANCER_TYPE
}

# Invalid NOFO reports
INVALID_NOFOS_REPORT = REPORTS_DIR +"/"+ "invalidNofoReport_" + QUALTRICS_TYPE + ".csv"
CORRECTED_INVALID_NOFOS_REPORT = REPORTS_DIR +"/"+ "invalidNofoReport_corrected.csv"
REVIEWED_NOFO_INPUT = INTERMED_DIR + QUALTRICS_VERSION +"/"+ "invalidNofoReport_reviewed.csv"

# Invalid Award reports 
INVALID_AWARD_REPORT = REPORTS_DIR +"/"+ "invalidAwardReport_" + QUALTRICS_TYPE + ".csv"
CORRECTED_INVALID_AWARD_REPORT = REPORTS_DIR +"/"+ "invalidAwardReport_corrected.csv"
REVIEWED_AWARD_INPUT = INTERMED_DIR + QUALTRICS_VERSION +"/"+ "invalidAwardReport_reviewed.csv"



# ---
# GRANTS CLEANING CONFIGURATION

# Earliest fiscal year (int) to use when gathering grants
API_EARLIEST_FISCAL_YEAR = 2000

# Grant/project fields to keep from NIH RePORTER results
API_FIELDS = [
    'project_num',
    'core_project_num',
    'appl_id',
    'fiscal_year',
    'project_title',
    'abstract_text',
    'pref_terms',
    'organization',
    'principal_investigators',
    'program_officers',
    'award_amount',
    'agency_ic_fundings',
    'award_notice_date',
    'project_start_date',
    'project_end_date',
    'opportunity_number', # Replaces full_foa
    'api_source_search' # Not from RePORTER. Added during data gathering step.
]

# List of name fields received from the API that require reformatting
API_NAME_FIELDS = [
    'principal_investigators',
    'program_officers',
]

# Field name containing nested JSON with NCI funding
API_AGENCY_FUNDING_FIELD = 'agency_ic_fundings'

# List of organization fields nested within organization field from API
API_ORG_FIELD = 'organization'
API_ORG_SUBFIELDS =   [
    'org_name',
    'org_city',
    'org_state',
    'org_country'
]

# Field name for abstract text
ABSTRACT_TEXT_FIELD = 'abstract_text'

# Dictionary of old:new column names. Rename API fields to match INS terms
# Any terms not included will remain as retrieved from API
API_FIELD_RENAMER = {
    "project_num": "grant_id", # also used as GRANT_ID_FIELDNAME
    "core_project_num": "queried_project_id",
    "appl_id": "application_id",
    "pref_terms": "keywords",
    "agency_ic_fundings": "nci_funded_amount"
}

# Define column for grant ID sorting
GRANT_ID_FIELDNAME = 'grant_id'

# Define name for new program ID field
PROGRAM_ID_FIELDNAME = 'program.program_id'



# ---
# SUMMARY STATISTICS CONFIGURATION

# Dict of grants fields of interest and how to aggregate each
STAT_AGG_FUNCS_BY_COL = {
    'api_source_search': 'nunique',
    'queried_project_id': 'nunique',
    'grant_id': 'nunique',
    'fiscal_year': 'min',
}
STAT_FISCALYEAR_COL = 'fiscal_year'
STAT_CORE_PROJECT_COL = 'queried_project_id'

# Summary statistic export filenames
STAT_GRANTS_BY_PROGRAM_FILENAME = REPORTS_GATHERED_DIR +"/"+ "grantsStatsByProgram.csv"
STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME = REPORTS_GATHERED_DIR +"/"+ "sharedProjectsByProgramPair.csv"



#---
# PROJECTS CONFIGURATION

# Report of shared project value validation
MISMATCHED_PROJECT_VALUES_REPORT = REPORTS_GATHERED_DIR +"/"+ "mismatchedProjectValuesReport.csv"



# ---
# PUBLICATIONS CONFIGURATION

# ICite bulk download csv.zip location
ICITE_FILENAME = INPUT_DIR +"icite/"+ ICITE_VERSION +"/"+ "icite_metadata.zip"

# Versioned directories for intermediates and outputs
TEMP_PUBLICATION_DIR = GATHERED_DIR +"/"+ "temp_pubmed_chunkfiles"
REMOVED_PUBLICATIONS = REPORTS_GATHERED_DIR +"/"+ "removedPublicationsReport.csv"
PROJECT_PMIDS = GATHERED_DIR +"/"+ "projectPMIDs.csv"
ICITE_PMID_DATA = GATHERED_DIR +"/"+ "icitePMIDData.csv"
MERGED_PMID_DATA = GATHERED_DIR +"/"+ "mergedPMIDData.csv"

# Publications output filepath
PUBLICATIONS_INTERMED_PATH = GATHERED_DIR +"/"+ "publication.csv"
PUBLICATIONS_OUTPUT_PATH = OUTPUT_GATHERED_DIR +"/"+ "publication.tsv"

# Earliest Publication year
PUBLICATION_YEAR_CUTOFF = 2000

# Temporary PubMed file chunksize
PUB_DATA_CHUNK_SIZE = 2000

# iCite columns of interest
ICITE_COLUMNS_TO_PULL = ['pmid','title','authors','year',
                         'citation_count','relative_citation_ratio']

# List of programs to exclude from downstream publication gathering
PROGRAMS_EXCLUDE_FROM_PUBS = ['ccdi']



# ---
# DATA PACKAGING CONFIGURATION

# Report subfolder for data packaing reports
PACKAGING_REPORTS = REPORTS_GATHERED_DIR +'/'+ 'packagingReports/'
REMOVED_DUPLICATES = PACKAGING_REPORTS + 'duplicate_' # Add datatype.csv in code
REMOVED_EARLY_PUBLICATIONS = PACKAGING_REPORTS + 'removedEarlyPublications.csv'

# Allowable difference between publication date and later project start date
PUB_PROJECT_DAY_DIFF = 365

# Generation of enum value lists
ENUM_PROGRAM_COLS = ['focus_area', 'cancer_type']
ENUM_PROGRAM_PATH = PACKAGING_REPORTS + 'program_enums.txt'

# Dictionary of columns and types to use in final data packaging
COLUMN_CONFIGS = {
    # Data type
    'program': {
        # Identifying node_id column name
        'node_id': 'program_id',
        # Idenfiying column name for relationship link
        'link_id': None,
        # Dict of old:new column names. Includes only columns to include in output
        'keep_and_rename': {
            'type': 'type',
            'program_id': 'program_id',
            'program_name': 'program_name',
            'program_acronym': 'program_acronym',
            'focus_area': 'focus_area',
            'cancer_type': 'cancer_type',
            'doc': 'doc',
            'contact_pi': 'contact_pi',
            'contact_pi_email': 'contact_pi_email',
            'contact_nih': 'contact_nih',
            'contact_nih_email': 'contact_nih_email',
            'nofo': 'nofo',
            'award': 'award',
            'program_link': 'program_link',
            'data_link': 'data_link',
        },
        # List of any list-like columns that need semicolon separators
        'list_like_cols': ['focus_area', 'cancer_type', 'doc'],
    },
    'grant': {
        'node_id': 'grant_id',
        'link_id': 'project.project_id',
        'keep_and_rename': {
            'type': 'type',
            'grant_id': 'grant_id',
            'queried_project_id': 'project.project_id',
            'application_id': 'application_id',
            'fiscal_year': 'fiscal_year',
            'project_title': 'grant_title',
            'abstract_text': 'abstract_text',
            'keywords': 'keywords',
            'org_name': 'org_name',
            'org_city': 'org_city',
            'org_state': 'org_state',
            'org_country': 'org_country',
            'principal_investigators': 'principal_investigators',
            'program_officers': 'program_officers',
            'award_amount': 'award_amount',
            'nci_funded_amount': 'nci_funded_amount',
            'award_notice_date': 'award_notice_date',
            'project_start_date': 'project_start_date',
            'project_end_date': 'project_end_date',
            'opportunity_number': 'opportunity_number', 
        },
        'list_like_cols': ['keywords', 'principal_investigators'],
        'datetime_cols': ['award_notice_date', 'project_start_date', 'project_end_date']
    },
    'project': {
        'node_id': 'project_id',
        'link_id': 'program.program_id',
        'keep_and_rename': {
            'type': 'type',
            'project_id': 'project_id',
            'program.program_id': 'program.program_id',
            'project_title': 'project_title',
            'abstract_text': 'abstract_text',
            'org_name': 'org_name',
            'org_city': 'org_city',
            'org_state': 'org_state',
            'org_country': 'org_country',
            'project_start_date': 'project_start_date',
            'project_end_date': 'project_end_date',
            'opportunity_number': 'opportunity_number',
        },
        'list_like_cols': ['opportunity_number'],
        'datetime_cols': ['project_start_date', 'project_end_date']
    },
    'publication': {
        'node_id': 'pmid',
        'link_id': 'project.project_id', 
        'keep_and_rename': {
            'type': 'type',
            'pmid': 'pmid',
            'coreproject': 'project.project_id',
            'title': 'title',
            'authors': 'authors',
            'publication_date': 'publication_date',
            'citation_count': 'cited_by',
            'relative_citation_ratio': 'relative_citation_ratio'
        },
        'list_like_cols': ['authors'],
        'html_tag_cols': ['title']
    }
}



# ---
# DATASETS CONFIGURATION

# dbGaP input file - CSV download of dbGaP search results
DBGAP_INPUT_CSV = INPUT_DIR + "dbgap/" + "study_" + DBGAP_CSV_VERSION + ".csv"

# dbGaP intermediate storage directory
DBGAP_INTERMED_DIR = INTERMED_DIR + "dbgap/" + DBGAP_CSV_VERSION + "/"

DBGAP_META_INTERMED_PATH = DBGAP_INTERMED_DIR + "dbgap_study_metadata.json"
DBGAP_SSTR_INTERMED_PATH = DBGAP_INTERMED_DIR + "dbgap_sstr.json"

DBGAP_PROCESSED_PATH = DBGAP_INTERMED_DIR + "dbgap_datasets.csv"

# dbGaP reports and error logs
DBGAP_REPORTS_DIR = "reports/dbgap/" + DBGAP_CSV_VERSION + "/"
DBGAP_SSTR_ERRORS = DBGAP_REPORTS_DIR + "api_errors_sstr.csv"
DBGAP_META_ERRORS = DBGAP_REPORTS_DIR + "api_errors_study_metadata.csv"

# dbGaP GPA/DOC input files
DBGAP_GPA_LIST = INPUT_DIR + "dbgap/gpa_tables/" + "gpa_study_table.csv"
DBGAP_GPA_DOC_LUT = INPUT_DIR + "dbgap/gpa_tables/" + "gpa_doc_lookup_table.csv"