# import
## batteries
import os
import sys
import time
import argparse
import requests
import functools
import multiprocessing as mp
import xml.etree.ElementTree as ET


# functions
def get_desc():
    desc = 'Get the NCBI lineage for >=1 taxonomy ID'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:
    Get NCBI lineages for >=1 taxonomy ID (in parallel).
    The entrez API is used for querying.

    If you have many thousands of IDs, then consider using taxID2LinTbl instead.
    """
    if subparsers:
        parser = subparsers.add_parser('taxID2lin', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    io = parser.add_argument_group('File input/output')
    io.add_argument('taxID', metavar='taxID', type=str,
                     help='Input table containing the taxonomy IDs, use "STDIN" if from STDIN')
    io.add_argument('-c', '--column', default=1,
                     help='Column number containing the accessions (default: %(default)s)')
    io.add_argument('-s', '--sep', default='\t',
                     help='Column separator (default: %(default)s)')
    io.add_argument('-x', '--header', default=False, action='store_true',
                     help='Header in input table? (default: %(default)s)')
    io.add_argument('-o', '--outfile', default='STDOUT',
                     help='Output file name; "STDOUT" if to STDOUT (default: %(default)s)')

    lin = parser.add_argument_group('Lineage')
    lin.add_argument('-l', '--levels', type=int, default=8,
                     help='Number of taxonomic levels (default: %(default)s)')                     

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-t', '--tries', type=int, default=3,
                      help='Number of tries to make each request. (default: %(default)s)')
    misc.add_argument('-p', '--procs', type=int, default=1,
                      help='Number of processors to use. (default: %(default)s)')

    # running test args
    if test_args:
        args = parser.parse_args(test_args)
        return args


def main(args=None):
    # Input
    if args is None:
        args = parse_args()

    # getting taxonomy IDs
    taxIDs = get_taxIDs(args.taxID, col_idx=args.column, sep=args.sep, header=args.header)

    # getting lineages
    procs = int(args.procs)
    if procs < 2:
        procs = None
    pool = mp.Pool(processes = procs)
    func = functools.partial(query_ncbi_lineage, levels=args.levels, tries=args.tries)
    lineages = pool.map(func, taxIDs)

    # writing lineages
    write_lineages(lineages, args.outfile, args.levels)


def write_lineages(lineages, outfile, levels):
    if outfile == 'STDOUT':
        outF = sys.stdout
    else:
        outF = open(outfile, 'r')

    header = ['taxID'] + ['rank_{}'.format(x+1) for x in range(levels)]
    outF.write('\t'.join(header) + '\n')

    for lin in lineages:
        outF.write('\t'.join(lin) + '\n')

    outF.close()
        

def get_taxIDs(infile, col_idx=1, sep='\t', header=False):
    col_idx = int(col_idx) - 1

    if infile == 'STDIN':
        inF = sys.stdin
    else:
        inF = open(infile)
    
    taxIDs = {}
    for i,line in enumerate(inF):
        if i == 0 and header == True:
            continue
        ID = line.strip().split(sep)[col_idx]
        if ID == '':
            continue
        try:
            taxIDs[ID] += 1
        except KeyError:
            taxIDs[ID] = 1

    if infile != 'STDIN':
        inF.close()

    # duplicates
    msg = 'WARNING: taxID "{}" was duplicated\n'
    for x,v in taxIDs.items():
        if v > 1:
            sys.stderr.write(msg.format(x))
            
    # ret
    return taxIDs.keys()
    

def query_ncbi_lineage(taxon_id, levels=9, tries=3):
    """Obtain the NCBI lineage for a taxon ID
    
    Parameters
    ----------
    taxon_id : int
        The taxon ID of interest
    
    Returns
    -------
    list or None
        Each taxon name or None if unable to retreive the taxon details
    """
    url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    # Define our parameters to use in our query
    params = {'db': 'taxonomy',  # We want to query the taxonomy database
              'id': taxon_id}    # We're requesting detail on the taxon ID specifically
    
    for i in range(tries):
        # Make the request
        r = requests.get(url, params=params)
    
        # Bail if we received a bad status
        if r.status_code != 200:
            if i < tries-1:
                msg = 'WARNING: status code = {} for taxID {}. Retrying\n'
                sys.stderr.write(msg.format(r.status_code, taxon_id))
                time.sleep(3)
            else:
                msg = 'WARNING: status code = {} for taxID {}. Giving up\n'
                sys.stderr.write(msg.format(r.status_code, taxon_id))
                return [taxon_id] + ['unclassified'] * levels
        else:
            break
    
    # NCBI returns XML, so we need to parse the "content" of our request into a usable structure
    tree = ET.fromstring(r.content)
    
    # We are only interested in the Lineage key within the tree. 
    # There is a fair bit of other data present within the structure, however.
    lineage = next(tree.iter('Lineage'))
    
    # The lineage, if we did get it, is delimited by "; "
    # So let's go ahead and parse that into a friendlier structure
    if lineage is not None:
        # splitting lineages
        lin = [v.strip() for v in lineage.text.split(';')]
        # expanding if needed
        if len(lin) > levels:
            lin = lin[0:levels]
        elif len(lin) < levels:
            lin = lin + ['unclassified'] * (levels - len(lin))
        return [taxon_id] + lin
        
    else:
        return None
