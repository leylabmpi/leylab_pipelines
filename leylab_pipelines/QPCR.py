# -*- coding: utf-8 -*-
# import
## batteries
from __future__ import print_function
import os
import sys
import argparse
from itertools import product
import string
## 3rd party
import numpy as np
import pandas as pd
## package
from leylab_pipelines import Utils
from leylab_pipelines import Fluent


# functions
def parse_args(test_args=None):
    # desc
    desc = 'Create robot commands for qPCR setup'
    epi = """DESCRIPTION:
    Create a worklist file for the TECAN Fluent robot for qPCR setup.
    The input is an exported plate layout from the BioRad qPCR software.
    The file format should be Excel or csv.
    The following columns should also be added to the table:
    * Sample labware (labware containing sample DNA/RNA; eg., "96-Well[001]")
    * Sample location (numeric; minimum of 1)
    * Sample volume (numeric)
    * MM name (Name of master mix; only used if you have multiple master mixes)
    * MM volume (Volume of master mix in PCR rxn)
    * Water volume (Volume of water in PCR rxn)
    
    Notes:
    * Sample locations in plates numbered are column-wise. 
    * The setup file (input table) MUST have a header (capitalization doesn't matter)
    * All volumes are in ul.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('setup', metavar='SetupFile', type=str,
                         help='An Excel or CSV file with experimental setup')
    groupIO.add_argument('--prefix', type=str, default='TECAN_qPCR',
                         help='Output file name prefix')
    groupIO.add_argument('--format', type=str, default=None,
                        help='File format (Excel or CSV). If not provided, the format is determined from the file extension') 

    # ## concentration file
    # conc = parser.add_argument_group('Concentation file')
    # conc.add_argument('--format', type=str, default=None,
    #                     help='File format (excel or tab). If not provided, the format is determined from the file extension') 
    # conc.add_argument('--header', action='store_false',
    #                     help='Header in the file? (Default: true)')
    # conc.add_argument('--labware', type=int, default=1,
    #                     help='Column containing the sample labware IDs') 
    # conc.add_argument('--location', type=int, default=2,
    #                     help='Column containing the sample location numbers')
    # conc.add_argument('--conc', type=int, default=3,
    #                     help='Column containing the sample concentrations')
    # conc.add_argument('--rows', type=str, default='all',
    #                   help='Which rows (not including header) of the column file to use ("all"=all rows; "1-48"=rows 1-48)')

    ## source sample plate
    src = parser.add_argument_group('Source sample plate')
    src.add_argument('--srctype', type=str, default='96-well',
                      choices=['96-well','384-well'],
                      help='Source sample plate labware type')
    src.add_argument('--srclabware', type=str,
                      default='96-well:96-well [001],384-well:384 Well[001],Tube:1x24 tube runner',
                      help='Choices for the source labware name base on --srctype')

    ## destination plate
    dest = parser.add_argument_group('Destination plate')
    dest.add_argument('--desttype', type=str, default='384-well',
                      choices=['96-well','384-well'],
                      help='Destination plate labware type')
    dest.add_argument('--destlabware', type=str,
                      default='96-well:96-well [002],384-well:384 Well[002]',
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
    # destination start
    args.destype = args.desttype.lower()
    # source labware
    args.srclabware = {x.split(':')[0].lower():x.split(':')[1] for x in args.srclabware.split(',')} 
    # destination labware
    args.destlabware = {x.split(':')[0].lower():x.split(':')[1] for x in args.destlabware.split(',')} 


def load_setup(input_file, file_format=None, header=0):
    # format
    if file_format is None:
        if input_file.endswith('.txt'):
            file_format = 'csv'
        elif input_file.endswith('.csv'):
            file_format = 'txt'
        elif input_file.endswith('.xls') or input_file.endswith('.xlsx'):
            file_format = 'excel'
    else:
        file_format = file_format.lower()
        
    # load via pandas IO
    if file_format == 'csv':
        df = pd.read_csv(input_file, sep=';', header=header)
    elif file_format == 'tab':
        df = pd.read_csv(input_file, sep='\t', header=header)
    elif file_format == 'excel':
        xls = pd.ExcelFile(input_file)
        df = pd.read_excel(xls, header=header)
    else:
        raise ValueError('Setup file not in usable format')

    # standarizing column IDs
    df.columns = [x.lower() for x in df.columns]

    # assert & return
    assert df.shape[1] > 1, 'Input file is only 1 column; wrong delimiter used?'    
    return df


def check_df_setup(df_setup):
    """Assertions of df_conc object formatting
    """
    # checking for column IDs
    col_IDs = ('row', 'column', 'sample type', 
               'sample labware', 'sample location', 'sample volume',
               'mm name', 'mm volume', 'water volume')
    msg = 'Column "{}" not found (captilization invariant)'
    for x in col_IDs:
        if not x in df_setup.columns:
            raise ValueError(msg.format(x))

    # checking sample locations (>=1)
    msg = 'ERROR (SetupFile, line={}): location is < 1'
    for i,loc in enumerate(df_setup['sample location']):
        if loc < 1:
            print(msg.format(i), file=sys.stderr)
    
    # checking sample conc
    msg = 'ERROR (setupFile, line={}): volume is <= 0'
    for i,vol in enumerate(df_setup['sample volume']):
        if vol <= 0.0:
            print(msg.format(i), file=sys.stderr)
    for i,vol in enumerate(df_setup['mm volume']):
        if vol <= 0.0:
            print(msg.format(i), file=sys.stderr)
    for i,vol in enumerate(df_setup['water volume']):
        if vol <= 0.0:
            print(msg.format(i), file=sys.stderr)
    

def _edit_src_labware(labware_type, labware_index):
    msg = 'Labware type "{}" not recognized'
    if isinstance(labware_type, str):
        labware_type = labware_type.lower()
        try: 
            x = labware_index[labware_type]
        except KeyError:
            raise KeyError(msg.format(labware_type))
    elif np.isnan(labware_type):
        return None
    else:
        raise ValueError('Logic error for labware type "{}"'.format(labware_type))        
    return x


def edit_src_labware(df_setup, src_labware_index, src_type='96-well'):
    """Editing source labware for TECAN gwl
    Changing ['sample labware'] columne
    """
    try:
        src_labware = src_labware_index[src_type]
    except KeyError:
        msg = 'src labware type "{}" not recognized'
        raise KeyError(msg.format(src_type))
    
    # changing sample labware
    func = lambda x: _edit_src_labware(x['sample labware'], src_labware_index)
    df_setup['sample labware'] = df_setup.apply(func, 1)


def plate2robot_loc(row_val, col_val, plate_type='96-well'):
    """Changing positioning from row (letter) and column (number)
    to just numeric position (column-wise) on the plate,
    which is needed for TECAN robot positioning. 
    Using index for identifying well number in plate
    [args]
    row_val: string
    col_vol: string
    """    
    # index for converting row to numeric
    idx = string.ascii_uppercase
    idx = {x:i+1 for i,x in enumerate(idx)}    
    row_val = idx[row_val]

    # getting location on plate
    msg = 'Destination location "{}" is out of range'
    if plate_type == '96-well':
        loc = (col_val - 1) * 8 + row_val
        assert loc > 0 and loc <= 96, msg.format(loc)
    elif plate_type == '384-well':
        loc = (col_val - 1) * 16 + row_val
        assert loc > 0 and loc <= 384, msg.format(loc)
    else:
        msg = 'Labware type "{}" not recognized'
        raise ValueError(msg.format(plate_type))
    return loc


def add_dest(df_setup, dest_labware_index, dest_type='96-well'):
    """Setting destination locations for samples & reagents
    Adding to df_conc:
      [dest_labware, dest_location]
    """
    # setting destination labware
    try:
        df_setup['dest_labware'] = dest_labware_index[dest_type]
    except KeyError:
        msg = 'Destination labware type "{}" not supported'
        raise KeyError(msg.format(dest_type))
    
    # setting destination location based on plate layout 
    func = lambda x: plate2robot_loc(x['row'], x['column'], plate_type=dest_type)
    df_setup['dest_location'] = df_setup.apply(func, 1)

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

    # Load input table
    df_setup = load_setup(args.setup, 
                          file_format=args.format)
    check_df_setup(df_setup, args)

    # adding sample source to setup table
    edit_src_labware(df_setup, 
                     src_type=args.srctype, 
                     src_labware_index=args.srclabware)

    # adding sample/reagent destinations to setup table
    add_dest(df_setup,
             dest_type=args.desttype,
             dest_labware_index=args.destlabware)
        
    # # Reordering dest if plate type is 384-well
    if args.desttype == '384-well':
         df_set = reorder_384well(df_setup, 'dest_location')
    
    # # Writing out gwl file
    # gwl_file = args.prefix + '.gwl'
    # with open(gwl_file, 'w') as gwlFH:
    #     ## Dilutant
    #     pip_dilutant(df_conc, outFH=gwlFH, src_labware=args.dlabware)
    #     ## Sample
    #     pip_samples(df_conc, outFH=gwlFH)

    # # Writing out table
    # conc_file = args.prefix + '_conc.txt'
    # df_conc.round(1).to_csv(conc_file, sep='\t', index=False)

    # # Create windows-line breaks formatted versions
    # gwl_file_win = Utils.to_win(gwl_file)
    # conc_file_win = Utils.to_win(conc_file)

    # # end
    # return (gwl_file, gwl_file_win, conc_file, conc_file_win)


# main
if __name__ == '__main__':
    pass


