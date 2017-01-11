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
    * Sample concentration (numeric value; units=ng/ul)
    
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
                     help='Target dilution concentration (ng/ul)')
    dil.add_argument('--minvolume', type=float, default=2.0,
                     help='Minimum sample volume to use')
    dil.add_argument('--maxvolume', type=float, default=30.0,
                     help='Maximum sample volume to use')
    dil.add_argument('--mintotal', type=float, default=10.0,
                     help='Minimum post-dilution total volume')

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
    args.destype = args.desttype.lower()
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

def dilution_volumes(df_conc, dilute_conc, min_vol, max_vol, 
                     min_total, dest_type='96-well'):
    """Setting the amoutn of sample to aliquot for dilution
    df_conc: pd.dataframe
    dilute_conc: concentration to dilute to 
    min_vol: min volume of sample to use
    max_vol: max total volume to use
    min_total: minimum total volume
    """
    # c1*v1 = c2*v2 
    # v2 = c1 * v1 / c2
    # v1 = c2 * v2 / c1

    # max well volume
    if dest_type == '96-well':
        max_well_vol = 280
    elif dest_type == '384-well':
        max_well_vol = 140
    else:
        raise ValueError('--desttype not recognized')

    # range of dilutions
    samp_vol_range = max_vol - min_vol
    target_total_vol = round(samp_vol_range / 2 + min_vol)
    
    # bare min for sample volume used
    df_conc['total_volume'] = [x * min_vol / dilute_conc for x in df_conc['conc']]
    if max(df_conc['total_volume']) > max_well_vol:
        msg = 'ERROR: dilution volume exceeds max possible well volume.'
        msg += ' Lower --minvolume or chane destination labware type.'
        raise ValueError(msg)
    
    # raising total_vols if low volume (if small dilution factor)
    df_conc.loc[df_conc.total_volume < min_total, 'total_volume'] = min_total
    # setting volumes
    ## v1 = c2*v2/c1
    ## sample_volume = dilute_conc * total_volume / conc 
    f = lambda row: dilute_conc * row['total_volume'] / row['conc']
    df_conc['sample_volume'] = df_conc.apply(f, axis=1)
    # dilutatant volume = total_volume - sample_volume
    f = lambda row: row['total_volume'] - row['sample_volume']
    df_conc['dilutant_volume'] = df_conc.apply(f, axis=1)
    
    # return
    return df_conc
        

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


def pip_dilutant(df_conc):
    """Writing worklist commands for aliquoting dilutant.
    Using 1-asp-multi-disp with 200 ul tips.
    Method:
    * calc max multi-dispense for 50 or 200 ul tips 
    """
    # making multi-disp object
    print('C;MasterMix')
    MD = Fluent.multi_disp()
    MD.SrcRackLabel = 'trough[001]'                        # user defined?
    MD.SrcPosition = 1                                     # need to set for all channels?
    MD.DestRackLabel = df_conc.dest_labware
    MD.DestPositions = df_conc.dest_location
    MD.Volume = df_conc.dilutant_volume                # NEED variable volumes!
    MD.NoOfMultiDisp = int(np.floor(180 / mmvolume))  # using 200 ul tips
    # writing
    print(MD.cmd() + '\n')


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

    # Determining dilution volumes
    df_conc = dilution_volumes(df_conc, 
                               dilute_conc=args.dilution,
                               min_vol=args.minvolume,
                               max_vol=args.maxvolume,
                               min_total=args.mintotal,
                               dest_type=args.desttype)
    
    # Adding destination data
    df_conc = add_dest(df_conc, 
                       dest_labware=args.destlabware,
                       dest_type=args.desttype,
                       dest_start=args.deststart)
    
    # Reordering dest if plate type is 384-well
    if args.desttype == '384-well':
        df_conc = reorder_384well(df_conc, 'dest_location')
    
    # Writing out gwl file
    

# main
if __name__ == '__main__':
    pass


