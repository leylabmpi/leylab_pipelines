# import
## batteries
import os
import sys
import glob
import tempfile
import argparse
import logging
## 3rd party
import dask.dataframe as dd
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
    Join tables with the dask python package, which enables parallel processing
    and fast execution.
    
    JOIN:
      To designate columns to join on, the --join input should be in the format:
      'X=X,Y=Y,...'
      For example, to join on 'accession' in table 1 with 'Acc' in table 2:
         --join 'accession=Acc'
      This can be extended to multiple column joins with a comma-sep list.

    OUTPUT:
      The table partitions will be written to separate temporary files, then
      joined into the final output file (designated with --outfile)
    """
    if subparsers:
        parser = subparsers.add_parser('join', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    io = parser.add_argument_group('Input/Output')
    io.add_argument('table', metavar='table', type=str, nargs=2,
                     help='File names of the tables to join on')
    io.add_argument('-j', '--join', default='1=1',
                     help='Columns to join on. See DESCRIPTION (default: %(default)s)')
    io.add_argument('-H', '--how', default='inner',
                     help='How to join the tables (default: %(default)s)')                    
    io.add_argument('-o', '--outfile', default='-',
                     help='Output file name; "-" if to STDOUT (default: %(default)s)')
    io.add_argument('-s', '--sep', default='\t',
                     help='Column separator (default: %(default)s)')
#    io.add_argument('-x', '--no-header', default=False, action='store_true',
#                     help='No header in input table (default: %(default)s)')

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-p', '--procs', default=1,
                     help='Number of processors to use. (default: %(default)s)')

    # running test args
    if test_args:
        args = parser.parse_args(test_args)
        return args
        

def parse_join(join_str):
    """Parsing joing string in form of 'X=X,Y=Y'
    """
    join_on = {}
    for x in join_str.split(','):
        y = x.split('=')
        if len(y) != 2:
            msg = '--join should be in format "X=X,Y=Y"'
        try:
            join_on['left'].append(y[0])
        except KeyError:
            join_on['left'] = [y[0]]
        try:
            join_on['right'].append(y[1])
        except KeyError:
            join_on['right'] = [y[1]]

    return(join_on)


def get_table(infile, sep='\t'):
    if infile.endswith('.gz'):
        compression = 'gzip'
    else:
        compression = None
    df = dd.read_csv(infile, sep=sep, compression=compression)
    return(df)


def write_tmp(df):
    # writing out temporary files
    tmpDir = tempfile.gettempdir()
    tempFile = next(tempfile._get_candidate_names())
    tempFile = os.path.join(tmpDir, tempFile + '*.txt')
    logging.info('writing temporary files to: {}'.format(tempFile))
    df.to_csv(tempFile, index=False, sep='\t')
    return tempFile


def write_table(infile, outfile, conn='w', header=True):
    if outfile == '-':
        with open(infile, 'r') as inF:
            if header == False:
                line = inF.readline()
            for line in inF:
                sys.stdout.write(line)
    else:
        with open(infile, 'r') as inF, open(outfile, conn) as outF:
            if header == False:
                line = inF.readline()
            for line in inF:
                outF.write(line)

def main(args=None):
    # Input
    if args is None:
        args = parse_args()
    args.procs = int(args.procs)

    # parsing the join arg
    join_on = parse_join(args.join)

    # creating table objects
    df1 = get_table(args.table[0], sep=args.sep)
    df2 = get_table(args.table[1], sep=args.sep)

    # joining (merging)
    df = dd.merge(df1, df2, left_on=join_on['left'], right_on=join_on['right'], how=args.how) #, npartitions=args.procs)
    
    # writing out temporary files
    tempFile_str = write_tmp(df)
    
    # cat files
    for i,infile in enumerate(glob.glob(tempFile_str)):
        if i == 0:
            conn = 'w'
            header = True
        else:
            conn = 'a'
            header = False        
        write_table(infile, args.outfile, conn=conn, header=header)
            
        
