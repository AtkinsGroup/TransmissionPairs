"""Functions for extracting data from the Los Alamos HIV Sequence Database.

Note on URLs:
    Los Alamos lookup for accession id, e.g. AB097177
    https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=AB097177

    NCBI accession id lookup, e.g. AB097177
    http://www.ncbi.nlm.nih.gov/nuccore/AB097177?report=graph&log$=seqview

    SSAM_SE_id lookup:
    https://www.hiv.lanl.gov/cgi-bin/BASIC_BLAST/basic_blast_pg.cgi?SSAM_SE_id=149746

    Patient ID lookup:
    https://www.hiv.lanl.gov/components/sequence/HIV/search/patient.comp?pat_id=9008
"""

import re
from collections import OrderedDict
import logging
import warnings

import bs4
import requests
import pandas as pd

__author__ = "Stephen Gaffney"
__copyright__ = "Stephen Gaffney"
__license__ = "gpl3"

_logger = logging.getLogger(__name__)


def load_cluster(cluster_id):
    """Obtains cluster name, description, patient names+ids, and accession names+ids.

    The IDs are required for URLs.

    Args:
        cluster_id (int): a cluster ID (as in URL id)

    Returns:
        data (dict): Dictionary of 'cluster name', 'description', 'patients', 'accessions'.
    """
    cluster_name_str = 'Cluster Name'
    patient_url_format = r"patient.comp\?pat_id=(\d+)"  # includes backslash escape
    genbank_url_format = r"query_one.comp\?se_id=(\d+)"

    col_dict = {'cluster_str': 'Cluster Name',
                'desc': 'Cluster Description',
                'patients': 'Patient(s)',
                'accessions': 'Accession(s)',
                }

    cluster_path = 'https://www.hiv.lanl.gov/components/sequence/HIV/search/cluster.comp?clu_id={}'
    url = cluster_path.format(cluster_id)
    soup = _get_soup_from_url(url)
    tables = soup('table')

    # Get main table
    main_table = [i for i in tables if cluster_name_str in i.text]
    assert(len(main_table) == 1), 'Failed to find main table: multiple have cluster name field'
    main_table = main_table[0]

    df = pd.read_html(str(main_table))[0]
    df = df.set_index([0]).transpose()
    cluster_str = df[col_dict['cluster_str']].iloc[0]
    desc = df[col_dict['desc']].iloc[0].strip(""" '" """)

    # GET PATIENT IDS
    patients = OrderedDict()
    patient_links = soup.find_all(href=re.compile(patient_url_format))
    for i in patient_links:
        patient_name = i.text.strip()
        # href is of form 'patient.comp?pat_id=38516'
        href = i['href']
        patient_id = int(re.findall(patient_url_format, href)[0])
        patients[patient_name] = patient_id

    # GET SEQUENCE IDS
    seqs = OrderedDict()
    seq_links = soup.find_all(href=re.compile(genbank_url_format))
    for i in seq_links:
        seq_name = i.text.strip()
        # href is of form 'patient.comp?pat_id=38516'
        href = i['href']
        seq_id = int(re.findall(genbank_url_format, href)[0])
        seqs[seq_name] = seq_id

    data = {
        'cluster_name': cluster_str,
        'description': desc,
        'patients': patients,
        'accessions': seqs
    }

    return data


def search_db(max_rec=100, virus='HIV-1', subtype='A1*', region='GENOME'):
    """Builds dataframe of Sequence DB records for given field selections.

    Returns:
        df (pandas.DataFrame): DataFrame with the following columns:
            row_id, blast, patient_id, accession, seq_name, subtype, country,
            sampling_year, genomic_region, seq_length, organism
    """
    test_form = {
                 'slave': subtype,
                 'Genomic Region': region,
                 'max_rec': max_rec,
                 'show_sql': 'on',
                 'LENGTH': '100',
                 'master': virus,
                 'submit': 'Search',
                 'action': 'search',
                 }
    main_cols = [
        'row_id',
        'blast',
        'patient_comb',
        'accession',
        'seq_name',
        'subtype',
        'country',
        'sampling_year',
        'genomic_region',
        'seq_length',
        'organism',
    ]
    url = 'https://www.hiv.lanl.gov/components/sequence/HIV/search/search.comp'
    soup = _get_soup_from_url(url, data=test_form)
    df = _get_df_from_soup(soup, col_headers=main_cols)
    return df


class Cluster:
    """Holds patients."""
    def __init__(self, cluster_id: int):
        self.cluster_id = cluster_id
        data = load_cluster(cluster_id)
        self.cluster_name = data['cluster_name']
        self.description = data['description']
        self.comb_patients = data['patients']  # not needed?
        self.comb_accessions = data['accessions']
        self.patient_dict = OrderedDict()
        for patient_code, patient_id in self.comb_patients.items():
            patient = Patient(patient_id, patient_code)
            self.patient_dict[patient_id] = patient

        desc_list = []
        acc_list = []
        for patient_id, patient in self.patient_dict.items():
            desc_list.append(patient.desc)
            # acc_df = pd.DataFrame.from_records(patient.accession_list, columns=['accession_name', 'accession_id'])
            # acc_df.insert(0, 'patient_id', patient_id)
            acc_list.append(patient.accession_df)
        desc_df = pd.concat(desc_list, axis=1)
        # desc_df.index.name = None
        acc_df = pd.concat(acc_list, axis=0, ignore_index=True)
        assert (acc_df.accession_id.value_counts().max() <= 1), "Duplicate accession issue."
        self.desc_df = desc_df
        self.acc_df = acc_df


class Patient:
    def __init__(self, patient_id, patient_code=None):
        """Holds various patient attributes including accessions with timepoint information."""
        self.patient_id = patient_id
        self.patient_code = patient_code
        data = extract_patient_info(patient_id)
        self.desc = data['desc']
        self.accession_list = data['accessions']
        self.clusters = data['clusters']
        self.accession_df = extract_patient_accession_timepoints(patient_id)
        if len(self.accession_list) != len(self.accession_df):
            values_df = set(self.accession_df.accession_id.values)
            values_list = set([i[0] for i in self.accession_list])
            values_df_only = values_df.difference(values_list)
            values_list_only = values_list.difference(values_df)
            diff_list = []
            if values_df_only:
                diff_list.append('Basic accession listing did not include {}.'.format(','.join(values_df_only)))
            if values_list_only:
                diff_list.append('Extended accession listing did not include {}.'.format(','.join(values_list_only)))
            diff_str = ' '.join(diff_list)
            _logger.error(("Basic accession listing for patient {p} contains {nb} accessions, while timepoints table "
                           "lists {ne} accessions. {diffs}")
                          .format(p=patient_id, nb=len(self.accession_list), ne=len(self.accession_df),
                                  diffs=diff_str))


def extract_patient_info(patient_id: int):
    """Get patient info dictionary from patient_id.

    Returns:
        data (dict): Dictionary with keys: desc, accessions, clusters
    """
    patient_info_url = "https://www.hiv.lanl.gov/components/sequence/HIV/search/patient.comp?pat_id={}".format(patient_id)
    soup = _get_soup_from_url(patient_info_url)
    ptables = pd.read_html(str(soup))

    """tables:
        0: tools
        1: text descriptors
        2: accessions
    """
    desc = ptables[1].iloc[:, :2]
    desc = desc.rename(columns={0: 'var', 1: 'val'})
    desc = desc.loc[:desc.query("var == 'Accession(s)'").index.values[0]-1].set_index('var')['val']
    desc.name = patient_id

    acc_urls = soup(href=re.compile('asearch/query_one'))
    accessions = []
    for a in acc_urls:
        accession_id = a.text.strip()
        se_id = int(re.findall('se_id=(\d+)', a['href'])[0])
        accessions.append((accession_id, se_id))

    cluster_urls = soup(href=re.compile('cluster.comp'))
    clusters = []
    for a in cluster_urls:
        cluster_name = a.text.strip()
        clu_id = int(re.findall('clu_id=(\d+)', a['href'])[0])
        clusters.append((cluster_name, clu_id))

    return {'desc': desc,
            'accessions': accessions,
            'clusters': clusters,
            }


def _soup_pager(url):
    soup = _get_soup_from_url(url)
    yield soup
    while _has_next_page_not_final(soup) or _has_next_page_is_final(soup):
        page_id = _get_results_page_id(soup)
        soup = _get_soup_from_url(url, data={'action Next.x': 1, 'action Next.y': 1, 'id': page_id})
        yield soup


def extract_patient_accession_timepoints(patient_id: int):
    """Load large sequence accession table with timepoint information for patient_id."""
    time_url = ("https://www.hiv.lanl.gov/components/sequence/HIV/search/d_search.comp"
                "?ssam_pat_id={}&ssam_postfirstsample_days=*&ssam_poststarttreatment_days=*"
                "&ssam_postendtreatment_days=*&ssam_postseroconv_days=*&ssam_postinfect_days=*&ssam_fiebig=*") \
        .format(patient_id)

    main_cols = [
        'row_id',
        'blast',
        'patient_comb',
        'accession_id',
        'seq_name',
        'subtype',
        'country',
        'sampling_year',
        'days_from_first_sample',
        'fiebig_stage',
        'days_from_treatment_end',
        'days_from_treatment_start',
        'days_from_infection',
        'days_from_seroconversion',
        'genomic_region',
        'seq_length',
        'organism',
    ]

    df_list = []
    for ind, soup in enumerate(_soup_pager(time_url)):
        if ind:
            _logger.info("Loading page {} for patient {}".format(ind + 1, patient_id))
        df = _get_df_from_soup(soup, col_headers=main_cols)
        df_list.append(df)

    final_cols = df_list[0].columns
    df = pd.concat(df_list, axis=0, ignore_index=True)
    df = df[final_cols]
    df.drop('patient_comb', axis=1, inplace=True)
    return df


def _get_soup_from_url(url, data=None):
    """Get BeautifulSoup object from HIV Database url, using POST request if data supplied."""
    # @TODO: find a way to get around "certificate verify failed" error without verify=False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if data:
            request = requests.post(url, data=data, verify=False)
        else:
            request = requests.get(url, verify=False)
    content = request.content
    soup = bs4.BeautifulSoup(content, features="lxml", from_encoding='utf8')
    return soup


def _get_df_from_soup(soup, col_headers=None):
    """Extract table from html soup and verify no other tables present.

    Args:
        soup: BeautifulSoup object that includes table.
        col_headers (list): column headers, for renaming columns (REQUIRED).

    Returns:
        df (pd.DataFrame): table of parsed data
    """
    if not col_headers:
        raise Exception("Column headers must be specified for table.")
    links = soup(href=re.compile('patient.comp'))
    table = links[0].find_parents('table')[0]
    # Verify that there was only one such table
    temp = {id(t.find_parent('table')) for t in links}
    assert (len(temp) == 1), "Multiple table matches with accession links"
    # Read main table
    df = pd.read_html(str(table))[0]
    # allow for double colspan in 2nd column producing extra parsed column
    if len(col_headers) == len(df.columns) - 1:
        col_headers.insert(2, 'blast2')
    df.columns = col_headers
    # handle case where first row contains headers
    if ''.join([str(i) for i in df.iloc[0, :2].values]) == '#Select':
        df = df.iloc[1:].copy()
    # Split combined patient identifier column
    patient_codes, patient_ids = zip(*df.patient_comb.apply(_get_patient_ids))
    df.insert(2, 'patient_id', patient_ids)
    df.insert(3, 'patient_code', patient_codes)
    ncbi_links = table(href=re.compile('nuccore'))
    df['pos'], df['ncbi_url'] = zip(*[_process_ncbi_link(i) for i in ncbi_links])
    blast_urls = pd.Series([i['href'] for i in table.findAll('a', {'href': re.compile('blast')})])
    ssam_ids = blast_urls.apply(_get_ssam_se_id)
    df.insert(5, 'blast_ssam_se_id', ssam_ids)
    return df


def _has_next_page_not_final(soup):
    """Test for additional results pages, and """
    input_obj = soup.find('input', title="Next")
    if input_obj:
        return True
    else:
        return False


def _has_next_page_is_final(soup):
    next_obj = soup.find('input', title="Next")
    last_obj = soup.find('input', title="Last")
    if last_obj and not next_obj:
        return True
    else:
        return False


def _get_results_page_id(soup):
    return soup.find('input', attrs={'name': 'id'})['value']


def _process_ncbi_link(link):
    """
    Args:
        link (:obj:`bs4.element.Tag`): contains NCBI link.
    """
    title = link.find('img')['title']
    href = link.attrs['href']
    match = re.match('Start: (\d+)\s+Stop: (\d+). Link to NCBI sequence viewer', title)
    pos = ':'.join(match.groups()) if match else title
    return pos, href


def _get_patient_ids(patient_comb: str):
    """Split combined patient identifier into patient code and patient id.

    Info:
        patient_id is unique, and used for patient lookups.
        https://www.hiv.lanl.gov/components/sequence/HIV/search/patient.comp?pat_id=9008
    """
    r = re.match('([^\(]*)\(([^\)]*)\)', patient_comb.strip())
    patient_code, patient_id = r.groups()
    patient_id = int(patient_id.strip())
    return patient_code, patient_id


def _get_ssam_se_id(blast_url: str):
    """Get SSAM_SE_id from blast url.

    Example: https://www.hiv.lanl.gov/cgi-bin/BASIC_BLAST/basic_blast_pg.cgi?SSAM_SE_id=149746
    """
    val = re.findall('\_SE\_id\=(\d+)', blast_url)[0]
    return val

