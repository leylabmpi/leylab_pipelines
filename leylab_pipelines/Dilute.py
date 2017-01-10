# -*- coding: utf-8 -*-
# import
## batteries
from __future__ import print_function
import os
import sys
import argparse
from itertools import product
## 3rd party
import numpy as np
import pandas as pd
## package
from leylab_pipelines import Utils
from leylab_pipelines import Fluent


# functions
def parse_args(test_args=None):
    # desc
    desc = 'Convert a mapping file to a NGS amplicon worklist file for the TECAN robot'
    epi = """DESCRIPTION:
    Create a worklist file for the TECAN Fluent robot for diluting samples.
    The input is an Excel or tab-delimited file with:
    * Sample labware  (eg., "96-Well[001]")
    * Sample location (numeric value; minimum of 1)
    * Sample concentration (numeric value)
    
    You can designate the columns for each value (see options).

    Sample locations in plates numbered are column-wise. 

    All volumes are in ul.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('concfile', metavar='ConcFile', type=str,
                         help='An excel or tab-delim file of concentrations')

    ## concentration file
    conc = parser.add_argument_group('Concentation file')
    conc.add_argument('--format', type=str, default=None,
                        help='File format (excel or tab). If not provided, the format is determined from the file extension') 
    conc.add_argument('--header', action='store_false',
                        help='Header in the file? (Default: true)')
    conc.add_argument('--labware', type=int, default=1,
                        help='Column containing the sample labware IDs') 
    conc.add_argument('--location', type=int, default=2,
                        help='Column containing the sample location numbers')
    conc.add_argument('--conc', type=int, default=3,
                        help='Column containing the sample concentrations')
    conc.add_argument('--rows', type=str, default='all',
                      help='Which rows (not including header) of the column file to use ("all"=all rows; "1-48"=rows 1-48)')

    ## dilution
    dil = parser.add_argument_group('Dilution')
    dil.add_argument('--dilution', type=float, default=1.0,
                        help='Target dilution concentration')
    dil.add_argument('--minvolume', type=float, default=1.0,
                        help='Minimum volume of sample to ')
    dil.add_argument('--maxvolume', type=float, default=100.0,
                        help='Maximum total final volume (sample + dilutant)')

    ## destination plate
    dest = parser.add_argument_group('Destination plate')
    dest.add_argument('--desttype', type=str, default='96-well',
                      choices=['96-well','384-well'],
                      help='Destination plate labware type')
    dest.add_argument('--deststart', type=int, default=1,
                      help='Start well number on destination plate')
    dest.add_argument('--destlabware', type=str,
                      default='96-well:96 Well[008],384-well:384 Well[004]',
                      help='Choices for the destination labware name base on --desttype')

    # parse & return
    if test_args:
        args = parser.parse_args(test_args)
    else:
        args = parser.parse_args()
    return args


def check_args(args):
    """Checking user input
    """
    # input table column IDs
    assert args.labware >= 1, '--labware must be >= 1'
    assert args.location >= 1, '--location must be >= 1'
    assert args.conc >= 1, '--conc must be >= 1'
    # input table row select
    args.rows = Utils.make_range(args.rows, set_zero_index=True)
    # dilution
    assert args.dilution >= 0.0, '--dilution must be >= 0'
    assert args.minvolume >= 0.0, '--minvolume must be >= 0'
    assert args.maxvolume > 0.0, '--maxvolume must be > 0'
    # destination start
    if args.desttype == '96-well':
        destlimit = 96
    elif args.desttype == '384-well':
        destlimit = 384
    if args.deststart < 1 or args.deststart > destlimit:
        msg = 'Destination start well # must be in range: 1-{}'
        raise ValueError(msg.format(destlimit))
    # destination labware
    args.destlabware = {x.split(':')[0]:x.split(':')[1] for x in args.destlabware.split(',')}    

def conc2df(concfile, row_select=None, file_format=None, header=True,
            labware_col=1, location_col=2, conc_col=3):
    """Loading a concentration file as a pandas dataframe
    """
    if header==True:
        header=0
    else:
        header=None
    # format
    if file_format is None:
        if concfile.endswith('.txt') or concfile.endswith('.csv'):
            file_format = 'tab'
        elif concfile.endswith('.xls') or concfile.endswith('.xlsx'):
            file_format = 'excel'
    else:
        file_format = file_format.lower()
        
    # load via pandas IO
    if file_format == 'tab':
        df = pd.read_csv(concfile, sep='\t', header=header)
    elif file_format == 'excel':
        xls = pd.ExcelFile(concfile)
        df = pd.read_excel(xls, header=header)
    else:
        raise ValueError('Concentration file not in usable format')

    # selecting particular rows
    if row_select is not None:
        df = df.iloc[row_select]

    # adding column info
    ## changing column names
    cols = df.columns.tolist()
    cols[labware_col-1] = 'labware'
    cols[location_col-1] = 'location'
    cols[conc_col-1] = 'conc'
    df.columns = cols
    
    # return
    return df

def check_df_conc(df_conc, args):
    """Assertions of df_conc object formatting
    """
    # checking sample locations (>=1)
    msg = 'ERROR (concfile, line={}): location is < 1'
    for i,loc in enumerate(df_conc['location']):
        if loc < 1:
            print(msg.format(i), file=sys.stderr)
    
    # checking sample conc
    msg = 'ERROR (concfile, line={}): concentration is <= 0'
    for i,sc in enumerate(df_conc['conc']):
        if sc <= 0.0:
            print(msg.format(i), file=sys.stderr)

def add_dest(df_conc, dest_labware_index, dest_type='96-well', dest_start=1):
    """Setting destination locations for samples & primers.
    Adding to df_conc:
      [dest_labware, dest_location]
    """
    dest_start= int(dest_start)
    try:
        dest_labware = dest_labware_index[dest_type]
    except KeyError:
        raise KeyError('Destination labware type not recognized')
    
    # adding columns
    df_conc['dest_labware'] = dest_labware
    dest_end = dest_start + df_conc.shape[0] 
    df_conc['dest_location'] = list(range(dest_start, dest_end))

    # return
    return df_conc

def reorder_384well(df, reorder_col):
    """Reorder values so that the odd, then the even locations are
    transferred. This is faster for a 384-well plate
    df: pandas.DataFrame
    reorder_col: column name to reorder
    """
    df['TECAN_sort_IS_EVEN'] = [x % 2 == 0 for x in df[reorder_col]]
    df.sort_values(by=['TECAN_sort_IS_EVEN', reorder_col], inplace=True)
    df = df.drop('TECAN_sort_IS_EVEN', 1)
    df.index = range(df.shape[0])
    return df


def main(args=None):
    # Input
    if args is None:
        args = parse_args()
    check_args(args)
    # Import
    df_conc = conc2df(args.concfile, 
                      file_format=args.format,
                      row_select=args.rows, 
                      header=args.header,
                      labware_col=args.labware,
                      location_col=args.location, 
                      conc_col=args.conc)
    check_df_conc(df_conc, args)
    # Adding destination data
    add_dest(df_conc, 
             dest_labware=args.destlabware,
             dest_type=args.desttype,
             dest_start=args.deststart)

    # Determining sample volume to aliquot

    # Determining dilutant volume to aliquot

    # Reordering dest if plate type is 384-well
    if args.desttype == '384-well':
        df_conc = reorder_384well(df_conc, 'dest_location')
    
    # Writing out gwl file


# main
if __name__ == '__main__':
    pass


