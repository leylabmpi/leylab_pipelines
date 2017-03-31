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
    parser.add_argument('-i', '--input', default=None, nargs='+',
                        help='>=1 NCBI taxonomy file name (default: %(default)s)') 
    parser.add_argument('-u', '--url', default='ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/',
                        help='Base url for downloading the NCBI taxonomy files (default: %(default)s)')
    parser.add_argument('-t', '--types', default=['gb','wgs'],
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



def get_db(base_url, DBs, outDir=None):
    """Getting NCBI DB file
    Saving to a temporary directory by default
    """
    logging.info('downloading NCBI taxonomy dump...')
    
    # urls for DB files to download
    DB_psbl = {'wgs' = 'nucl_wgs.accession2taxid.gz',
               'gb' = 'nucl_est.accession2taxid.gz',
               'est' = 'nucl_est.accession2taxid.gz',
               'gss' = 'nucl_gss.accession2taxid.gz')

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
    

#--- getting taxIDs ---#

# # determining taxonomic IDs
# accs = list(seqs['accession'])
# taxIDs = []
# with gzip.open(acc2taxID, 'rb') as inF:
#     for line in inF:
#         line = line.rstrip().split('\t')
#         if line[0] in accs:
#             taxIDs.append(line[2])
#         else:
#             taxIDs.append(None)

# taxIDs[0:3]

def main(args=None):
    # Input
    if args is None:
        args = parse_args()

    # data downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    if args.outdir is not None or args.files is None:
        print('test')
#        files = get_taxdump(args.url, args.outdir)
#        args.nodes = files['nodes']
#        args.names = files['names']
        
