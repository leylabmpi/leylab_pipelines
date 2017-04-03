# import
## batteries
import re
import os
import sys
import gzip
import tempfile
import argparse
import logging
import urllib
## 3rd party
import pandas as pd
## package
from leylab_pipelines import Utils 

# logging
logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

# functions
def get_desc():
    desc = 'Get NCBI taxonomy IDs from NCBI accessions'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:

    TO CONVERT ACCESSION TO TAX_ID: 
      see ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    """
    if subparsers:
        parser = subparsers.add_parser('acc2taxID', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    parser.add_argument('accessions', metavar='accessions', type=str,
                        help='Input table containing the acessions. Use "-" for STDIN')
    parser.add_argument('-c', '--column', default=1,
                        help='Column number containing the accessions (default: %(default)s)')
    parser.add_argument('-s', '--sep', default='\t',
                        help='Column separator (default: %(default)s)')
    parser.add_argument('-x', '--no-header', default=False, action='store_true',
                        help='No header in input table (default: %(default)s)')
    parser.add_argument('-t', '--tax', default=None, nargs='+',
                        help='>=1 NCBI taxonomy file name (default: %(default)s)') 
    parser.add_argument('-u', '--url', default='ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/',
                        help='Base url for downloading the NCBI taxonomy files (default: %(default)s)')
    parser.add_argument('-T', '--types', default=['gb','wgs'], nargs='+',
                        help='>=1 DBs to download if no input files provided (default: %(default)s)') 
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


def load_acc_df(infile, sep='\t', no_header=False):
    logging.info('loading accessions...')

    if no_header is True:
        header = None
    else:
        header = 0
    if infile == '-':
        df = pd.read_csv(sys.stdin, sep=sep, header=header)
    else:
        df = pd.read_csv(infile, sep=sep, header=header)
    return df


def make_db_url(base_url, db):
    # urls for DB files to download
    DB_psbl = {'wgs' : 'nucl_wgs.accession2taxid.gz',
               'gb' : 'nucl_est.accession2taxid.gz',
               'est' : 'nucl_est.accession2taxid.gz',
               'gss' : 'nucl_gss.accession2taxid.gz'}
    
    try:
        x = DB_psbl[db]
    except KeyError:
        raise KeyError('DB type "{}" not recognized'.format(db))
    url = base_url + '/' + x 
    return [url, x]


def download_file(url, fileName, outDir=None):
    if outDir is None:
        outDir = tempfile.gettempdir()
    dbFile = os.path.join(outDir, fileName)

    if sys.version_info[0] >= 3:
        dbFile, headers = urllib.request.urlretrieve(url, filename=dbFile)
    else:
        dbFile, headers = urllib.urlretrieve(url, filename=dbFile)        

    logging.info('downloaded file: {}'.format(dbFile))
    return dbFile

    

def get_tax_db(base_url, DBs, outDir=None):
    """Getting NCBI DB file
    Saving to a temporary directory by default
    """
    logging.info('downloading NCBI taxonomy dump...')
    
    # url(s)
    urls = [make_db_url(base_url, x) for x in DBs]

    # downloading
    dmpFiles = [download_file(url[0], url[1], outDir) for url in urls]
                
    ## checking for existence
    for F in dmpFiles:
        if not os.path.isfile(F):
            raise ValueError('Cannot find file: {}'.format(F))

    return dmpFiles
    

def acc_to_taxID(db_file, df_acc, column=1):
    # accs 
    column = int(column) - 1
    accs = set(df_acc.iloc[:,column])
    # determining taxonomic IDs
    if db_file.endswith('.gz'):
        inF = gzip.open(db_file, 'rb') 
    else:
        inF =  open(db_file, 'r') 
    taxIDs = np.empty(df_acc.shape[0])
    for i,line in enumerate(inF):
        line = line.rstrip().split('\t')
        if line[0] in accs:
            taxIDs[i] = line[2])
r# # add column to df


def main(args=None):
    # Input
    if args is None:
        args = parse_args()

    # load accessions
    df_acc = load_acc_df(args.accessions, args.sep, args.no_header)

    # data downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    if args.outdir is not None or args.tax is None:
        args.tax = get_tax_db(args.url, args.types, args.outdir)

    # getting taxIDs for accessions
    for db_file in args.tax:
        acc_to_taxID(db_file, df_acc, column=args.column)
