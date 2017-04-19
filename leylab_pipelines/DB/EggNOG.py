# import
## batteries
import os
import sys
import argparse
import requests
import warnings
import functools
import multiprocessing as mp


# functions
def get_desc():
    desc = 'Get EggNOG data with the REST API'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:
    Download EggNOG database data (in parallel) with the EggNOG API.

    To find EggNOG group names, you can use the EggNOG website (http://eggnogdb.embl.de/#/app/home)
    """
    if subparsers:
        parser = subparsers.add_parser('eggnog', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    io = parser.add_argument_group('Input/Output')
    io.add_argument('nog_name', metavar='nog_name', type=str, nargs='+',
                     help='>= 1 valid EggNOG group name or files(s) of names (i.e. ENOG410ZSWV or COG0575)')
    io.add_argument('-a', '--attribute', default='fasta', 
                    choices=['fasta', 'raw_alg', 'trimmed_alg', 'tree', 'go_terms', 'domains'],
                    help='attribute (default: %(default)s)')    
    io.add_argument('-o', '--outdir', default='.',
                    help='Output directory (default: %(default)s)')                        

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-p', '--procs', type=int, default=1,
                      help='Number of processors to use (default: %(default)s)')

    # running test args
    if test_args:
        args = parser.parse_args(test_args)
        return args


def main(args=None):
    # args
    if args is None:
        args = parse_args()
    args.outdir = os.path.abspath(args.outdir)

    # getting nog_names
    nog_names = get_nog_names(args.nog_name)
    nog_names = filter_nog_names(nog_names)

    # writing content
    procs = int(args.procs)
    if procs < 2:
        procs = None
    pool = mp.Pool(processes = procs)
    func = functools.partial(write_content, attribute=args.attribute, outdir=args.outdir)
    outfiles = pool.map(func, nog_names)

    # status
    for x in outfiles:
        msg = 'File written: {}\n'
        sys.stderr.write(msg.format(x))

    
def read_nog_names(infile):
    nog_names = []
    with open(infile) as inF:        
        msg = 'Reading input file: {}\n'
        sys.stderr.write(msg.format(infile))
        for line in inF:
            line = line.rstrip().split('\t')[0]
            if line == '':
                continue
            nog_names.append(line)
    return nog_names


def get_nog_names(nog_name):
    nog_names = []
    for x in nog_name:
        try:
            nog_names += read_nog_names(x)
        except IOError:
            nog_names.append(x)

    return nog_names

def filter_nog_names(nog_names):
    # filtering
    filt = []
    for x in nog_names:
        if not (x.startswith('COG') or x.startswith('ENOG')):
            msg = 'WARNING: nog_name "{}" not recognized. Skipping\n'
            sys.stderr.write(msg.format(x))
        else:
            filt.append(x)

    if len(filt) == 0:
        sys.stderr.write('No acceptible nog names. Quitting\n')
        sys.exit(1)
                        
    return filt
        
            

def get_ext(attribute):
    psbl_ext = {'fasta' : '.fa',
                'raw_alg' : '.fa',
                'trimmed_alg' : '.fa',
                'tree' : '.nwk',
                'go_terms' : '.txt',
                'domains' : '.txt'}
    try:
        return psbl_ext[attribute.lower()]
    except KeyError:
        return ''


def write_content(nog_name, attribute, outdir):
    # creating output file name
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile = '_'.join([nog_name, attribute]) + get_ext(attribute)
    outfile = os.path.join(outdir, outfile)
    
    # attribute
    if attribute == 'go_terms' or attribute == 'domains':
        dataformat = 'json'
    else:
        dataformat = 'file'

    # request
    ## url
    url = 'http://eggnogapi.embl.de/nog_data/{}/{}/{}'
    url = url.format(dataformat, attribute, nog_name)
    print(url)
    ## checking request
    r = requests.get(url)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    ## writing
    with open(outfile, 'w') as outF:
        outF.write(r.text)
    
    return outfile
