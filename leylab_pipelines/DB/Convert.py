# import
## batteries
import argparse
import requests
import sys
import re
import warnings
import multiprocessing as mp
try:
    from urllib import urlencode
except ImportError:
    from urllib.parse import urlencode
## 3rd party
import xmltodict


# functions
def get_desc():
    desc = 'Convert between various gene or protein IDs'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:
    Convert between various gene or protein IDs using entrez API requests.

    Entrez requests can be done in parallel. 

    IDs:
      If the IDs are provided in a table, select the column delimiter and 
      column number (1-indexed). 
      The list of IDs can be provided via STDIN by using `STDIN`
    """
    if subparsers:
        parser = subparsers.add_parser('convert', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    io = parser.add_argument_group('Input/Output')
    io.add_argument('IDs', metavar='IDs', type=str,
                     help='A file containing IDs. See description.')
    io.add_argument('-m', '--method', default='Accession2Taxonomy',
                    choices=['Accession2Taxonomy', 
                             'EntrezGene2Uniprot',
                             'Uniprot2EntrezGene',
                             'HGNC2EntrezGene',
                             'EntrezTranGene2EntrezGene'],
                    help='ID conversion method. (default: %(default)s)')
    io.add_argument('-s', '--sep', default='\t',
                    help='Column separator (default: %(default)s)')
    io.add_argument('-c', '--column', type=int, default=1,
                    help='Column containing the IDs (default: %(default)s)')
    io.add_argument('-o', '--outfile', default='-',
                     help='Output file name; "-" if to STDOUT (default: %(default)s)')
    io.add_argument('-e', '--email', default='dummyemail@dummybunny.info',
                     help='user email address (default: %(default)s)')

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-p', '--procs', type=int, default=1,
                      help='Number of processors to use (default: %(default)s)')

    # running test args
    if test_args:
        args = parser.parse_args(test_args)
        return args


def main(args=None):
    # Input
    if args is None:
        args = parse_args()

    # reading IDs
    IDs = read_IDs(args.IDs, sep=args.sep, column=args.column)

    # conversion    
    IDs = get_conversion(IDs, args.email, args.method, procs=args.procs)
    
    # writing IDs
    write_IDs(IDs, args.outfile)


def write_IDs(IDs, outfile='-'):
    """writing out IDs
    """
    if outfile == '-':
        outF = sys.stdout
    else:
        outF = open(outfile, 'w')

    header = '\t'.join(['orig_ID', 'new_ID'])
    outF.write(header + '\n')
    for x in IDs:
        outF.write('\t'.join(x) + '\n')

    if outfile != '-':
        outF.close()    
    

def get_conversion(IDs, email, method, procs=1):
    conv = Conversion(email)

    if procs < 2:
        procs = None
    pool = mp.Pool(processes = procs)

    if method.lower() == 'accession2taxonomy':
        IDs = pool.map(conv.convert_accession_to_taxid, IDs)
    elif method.lower() == 'entrezgene2uniprot':
        IDs = pool.map(conv.convert_entrez_to_uniprot, IDs)
    elif method.lower() == 'uniprot2entrezgene':
        IDs = pool.map(convert_uniprot_to_entrez, IDs)
    elif method.lower() == 'hgnc2entrezgene':
        IDs = pool.map(convert_hgnc_to_entrez, IDs)
    elif method.lower() == 'entreztrangene2entrezgene':
        IDs = pool.map(convert_ensembl_to_entrez, IDs)
    else:
        msg = 'Method "{}" not recognized'
        raise IOError(msg.format(method))
    return IDs
    

def read_IDs(infile, sep='\t', column=1):
    """Reading in IDs
    """
    column = column - 1
    if infile == 'STDIN':
        inF = sys.stdin
    else:
        inF = open(infile, 'r')

    IDs = []
    for line in inF:
        line = line.rstrip().split(sep)
        IDs.append(line[column])

    if infile != 'STDIN':
        inF.close()
    return IDs


class Conversion(object):
    def __init__(self, email):
        """email is required
        """
        self.params = {}
        self.email = email
        self.params['tool'] = 'PyEntrez'
        if re.match(r"[^@]+@[^@]+\.[^@]+", self.email):
            pass
        else:
            raise ValueError("Enter a valid Email Address")
        self.params["email"] = email
        self.options = urlencode(self.params, doseq=True)
        
    def convert_ensembl_to_entrez(self, ensembl):
        """Convert Ensembl Id to Entrez Gene Id
        """        
        if 'ENST' in ensembl:
                pass
        else:
            raise(IndexError)
        # Submit resquest to NCBI eutils/Gene database
        server = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + self.options + "&db=gene&term={0}".format(ensembl)
        r = requests.get(server, headers={"Content-Type": "text/xml"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        # Process Request
        response = r.text
        info = xmltodict.parse(response)
        try:
            geneId = info['eSearchResult']['IdList']['Id']
        except TypeError:
            sys.stderr.write('WARNING: No Entrez ID for "{}"\n'.format(ensembl))
            geneId = None
        return [ensembl, geneId]

    def convert_hgnc_to_entrez(self, hgnc):
        """Convert HGNC Id to Entrez Gene Id
        """
        entrezdict = {}
        server = "http://rest.genenames.org/fetch/hgnc_id/{0}".format(hgnc)
        r = requests.get(server, headers={ "Content-Type" : "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        response = r.text
        info = xmltodict.parse(response)
        try:
            for data in info['response']['result']['doc']['str']:
                if data['@name'] == 'entrez_id':
                    entrezdict[data['@name']] = data['#text']
                if data['@name'] == 'symbol':
                    entrezdict[data['@name']] = data['#text']
        except KeyError:
            sys.stderr.write('WARNING: No Entrez ID for "{}"\n'.format(hgnc))            
            entrezdict = None
        try:
            entrezdict = entrezdict['entrez_id']
        except (KeyError, TypeError):
            entrezdict = None

        return [hgnc, entrezdict]

    def convert_entrez_to_uniprot(self, entrez):
        """Convert Entrez Id to Uniprot Id
        """
        server = "http://www.uniprot.org/uniprot/?query=%22GENEID+{0}%22&format=xml".format(entrez)
        r = requests.get(server, headers={ "Content-Type" : "text/xml"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        response = r.text
        info = xmltodict.parse(response)
        try:
            data = info['uniprot']['entry']['accession'][0]
        except TypeError:
            data = info['uniprot']['entry'][0]['accession'][0]
        return [entrez, data]

    def convert_uniprot_to_entrez(self, uniprot):
        """Convert Uniprot Id to Entrez Id
        """
        # Submit request to NCBI eutils/Gene Database
        server = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + self.options + "&db=gene&term={0}".format(uniprot)
        r = requests.get(server, headers={ "Content-Type" : "text/xml"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        # Process Request
        response = r.text
        info = xmltodict.parse(response)
        geneId = info['eSearchResult']['IdList']['Id']
        # check to see if more than one result is returned
        # if you have more than more result then check which Entrez Id returns the same uniprot Id entered.
        if len(geneId) > 1:
            for x in geneId:
                c = self.convert_entrez_to_uniprot(x)[1]
                c = c.lower()
                u = uniprot.lower()
                if c==u:
                    return [uniprot, x]
        else:
            return [uniprot, geneId]

    def convert_accession_to_taxid(self, accessionid):
        """Convert Accession Id to Tax Id
        """
        # Submit request to NCBI eutils/Taxonomy Database
        server = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + self.options + "&db=nuccore&id={0}&retmode=xml".format(accessionid)
        r = requests.get(server, headers={ "Content-Type" : "text/xml"})
        if not r.ok:
            try:
                r.raise_for_status()
            except requests.exceptions.HTTPError:
                sys.stderr.write('WARNING: No taxonomy ID for "{}"\n'.format(accessionid)) 
            return [accessionid, None]
        # Process Request
        response = r.text
        records = xmltodict.parse(response)
        try:
            for i in records['GBSet']['GBSeq']['GBSeq_feature-table']['GBFeature']['GBFeature_quals']['GBQualifier']:
                for key, value in i.items():
                    if value == 'db_xref':
                        taxid = i['GBQualifier_value']
                        taxid = taxid.split(':')[1]
                        return [accessionid, taxid]
        except:
            for i in records['GBSet']['GBSeq']['GBSeq_feature-table']['GBFeature'][0]['GBFeature_quals']['GBQualifier']:
                for key, value in i.items():
                    if value == 'db_xref':
                        taxid = i['GBQualifier_value']
                        taxid = taxid.split(':')[1]
                        return [accessionid, taxid]
        return [accessionid, None]
