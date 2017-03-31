# import
## batteries
import re
import os
import sys
import gzip
import tempfile
import multiprocessing
import argparse
import logging
import urllib
import tarfile
## 3rd party
import pandas as pd
## package
from leylab_pipelines import Utils 

# logging
logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

# functions
def get_desc():
    desc = 'Convert an NCBI taxonomy (taxdump) into lineages (taxID<-->lineage)'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:
    Convert the NCBI taxonomy dump to a lineage table.
    The lineage table can be used to map taxonomy IDs to taxonomy. 

    If the specific taxonomy dump files are not provided, this command
    will download the full dump file, uncompress it, and use the needed
    files. By default, the dump files are downloaded to a temporary directory.
    You can keep the taxonomy dump by using the --outdir option.

    TO CONVERT ACCESSION TO TAX_ID: 
      see ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    """
    if subparsers:
        parser = subparsers.add_parser('taxID2lin', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    parser.add_argument('--nodes', default=None,
                        help='nodes.dmp file path. (default: %(default)s)')
    parser.add_argument('--names', default=None,
                        help='names.dmp file path. (default: %(default)s)')
    parser.add_argument('-u', '--url', default='ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz',
                        help='URL for downloading taxonomy dump. (default: %(default)s)')
    parser.add_argument('-o', '--outfile', default='NCBI_taxID2lin.txt',
                        help='Output file for lineage table (default: %(default)s)')
    parser.add_argument('-d', '--outdir', default=None,
                        help='Output directory for the taxonomy dump download. (default: %(default)s)')
    parser.add_argument('-p', '--procs', default=1,
                        help='Number of processors to use. (default: %(default)s)')

    # running test args
    if test_args:
        args = parser.parse_args(test_args)
        return args


def strip(str_):
    '''
    :param str_: a string
    '''
    return str_.strip()


def load_nodes(nodes_file):
    '''
    load nodes.dmp and convert it into a pandas.DataFrame
    '''
    df = pd.read_csv(nodes_file, sep='|', header=None, index_col=False,
                     names=[
                         'tax_id',
                         'parent_tax_id',
                         'rank',
                         'embl_code',
                         'division_id',
                         'inherited_div_flag',
                         'genetic_code_id',
                         'inherited_GC__flag',
                         'mitochondrial_genetic_code_id',
                         'inherited_MGC_flag',
                         'GenBank_hidden_flag',
                         'hidden_subtree_root_flag',
                         'comments'
                     ])

    # To get rid of flanking tab characters
    df['rank'] = df['rank'].apply(strip)
    df['embl_code'] = df['embl_code'].apply(strip)
    df['comments'] = df['comments'].apply(strip)
    return df


def load_names(names_file):
    '''
    load names.dmp and convert it into a pandas.DataFrame
    '''
    df = pd.read_csv(names_file, sep='|', header=None, index_col=False,
                     names=[
                         'tax_id',
                         'name_txt',
                         'unique_name',
                         'name_class'
                     ])
    df['name_txt'] = df['name_txt'].apply(strip)
    df['unique_name'] = df['unique_name'].apply(strip)
    df['name_class'] = df['name_class'].apply(strip)

    sci_df = df[df['name_class'] == 'scientific name']
    sci_df.reset_index(drop=True, inplace=True)
    return sci_df


def to_dict(lineage):
    """
    convert the lineage into a list of tuples in the form of

    [
        (tax_id1, rank1, name_txt1),
        (tax_id2, rank2, name_txt2),
        ...
    ]

    to a dict
    """
    dd = {}
    num_re = re.compile('[0-9]+')
    len_lineage = len(lineage)
    for k, __ in enumerate(lineage):
        tax_id, rank, name_txt = __
        # use the last rank as the tax_id, whatever it is, genus or species.
        if k == len_lineage - 1:
            dd['tax_id'] = tax_id

        # e.g. there could be multiple 'no rank'
        numbered_rank = rank
        while numbered_rank in dd:
            # print __, numbered_rank
            search = num_re.search(numbered_rank)
            if search is None:
                count = 1
            else:
                count = int(search.group()) + 1
            numbered_rank = '{0}{1}'.format(rank, count)
        dd[numbered_rank] = name_txt
    return dd


def find_lineage(tax_id):
    if tax_id % 50000 == 0:
        logging.debug('working on tax_id: {0}'.format(tax_id))
    lineage = []
    while True:
        rec = TAXONOMY_DICT[tax_id]
        lineage.append((rec['tax_id'], rec['rank'], rec['name_txt']))
        tax_id = rec['parent_tax_id']

        if tax_id == 1:
            break

    # reverse results in lineage of Kingdom => species, this is helpful for
    # to_dict when there are multiple "no rank"s
    lineage.reverse()
    return to_dict(lineage)


def get_taxdump(url, outDir=None):
    """Getting taxdump file from NCBI. 
    Saving to a temporary directory. 

    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    """
    logging.info('downloading NCBI taxonomy dump...')
    
    # downloading
    if outDir:
        # saving taxdump
        dmpFile = os.path.join(outDir, 'taxdump.tar.gz')
    else:
        # taxdump written to temporary directory
        outDir = tempfile.gettempdir()
        dmpFile = None
    if sys.version_info[0] >= 3:
        dmpFile, headers = urllib.request.urlretrieve(url, filename=dmpFile)
    else:
        dmpFile, headers = urllib.urlretrieve(url, filename=dmpFile)        
            
    # uncompressing
    logging.info('uncompressing NCBI taxonomy dump file...')
    tar = tarfile.open(dmpFile, "r:gz")
    tar.extractall(path=outDir)
    tar.close()
    logging.info('files uncompressed to: {}'.format(outDir))

    # target file names
    files = {'nodes' : os.path.join(outDir, 'nodes.dmp'),
             'names' : os.path.join(outDir, 'names.dmp')}
    ## checking for existence
    for k,v in files.items():
        if not os.path.isfile:
            raise ValueError('Cannot find file: {}'.format(v))
    # ret
    return files
    

def main(args=None):
    # Input
    if args is None:
        args = parse_args()

#    print( int(args.procs).__class__ ); 
#    print( multiprocessing.cpu_count().__class__ );
#    exit()

#    ncpus = multiprocessing.cpu_count()
#    ncpus = int(args.procs)
#    logging.info('using {0} cpus to find lineages for all tax ids'.format(ncpus))
#    pool = multiprocessing.Pool(ncpus)
#    pool.close()
    
#    exit;
#    sys.exit;

    # data downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    if args.outdir is not None or args.nodes is None or args.names is None:
        files = get_taxdump(args.url, args.outdir)
        args.nodes = files['nodes']
        args.names = files['names']
        
    nodes_df = load_nodes(args.nodes)
    names_df = load_names(args.names)
    df = nodes_df.merge(names_df, on='tax_id')
    df = df[['tax_id', 'parent_tax_id', 'rank', 'name_txt']]
    df.reset_index(drop=True, inplace=True)
    logging.info('# of tax ids: {0}'.format(df.shape[0]))
    # log summary info about the dataframe
    df.info()

    # force to use global variable TAXONOMY_DICT because map doesn't allow
    # passing extra args easily
    global TAXONOMY_DICT
    logging.info('generating TAXONOMY_DICT...')
    TAXONOMY_DICT = dict(zip(df.tax_id.values, df.to_dict('records')))

#    ncpus = multiprocessing.cpu_count()
    ncpus = int(args.procs)
    logging.info('using {0} cpus to find lineages for all tax ids'.format(ncpus))
    pool = multiprocessing.Pool(ncpus)
    lineages_dd = pool.map(find_lineage, df.tax_id.values)
    pool.close()

    logging.info('generating a dictionary of lineages information...')
    dd_for_df = dict(zip(range(len(lineages_dd)), lineages_dd))

    logging.info('generating lineages_df...')
    lineages_df = pd.DataFrame.from_dict(dd_for_df, orient='index')
    lineages_df.sort_values('tax_id', inplace=True)

    logging.info('writing lineages to: {}'.format(args.outfile))
    cols = ['tax_id',
            'superkingdom',
            'phylum',
            'class',
            'order',
            'family',
            'genus',
            'species']
     
    other_cols = sorted([__ for __ in lineages_df.columns if __ not in cols])
    output_cols = cols + other_cols
    lineages_df.to_csv(args.outfile, index=False, columns=output_cols, sep='\t')

