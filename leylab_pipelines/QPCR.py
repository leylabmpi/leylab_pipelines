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
    The file format should be Excel or CSV.
    Just create a plate layout for your experimet, then export and add some needed info:
      The following columns should also be added to the table:
      * "Sample labware"  (labware containing sample DNA/RNA; eg., "96-Well[001]")
      * "Sample location"  (numeric; minimum of 1)
      * "Sample volume"  (numeric)
      * "MM name"  (Name of master mix; this allows for multiple mastermixes)
      * "MM volume"  (Volume of master mix in PCR rxn)
      * "Water volume"  (Volume of water in PCR rxn)
    
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
    msg = 'ERROR (setupFile, line={}): volume is < 0'
    for i,vol in enumerate(df_setup['sample volume']):
        assert np.isnan(vol) or vol >= 0.0, msg.format(i)
    for i,vol in enumerate(df_setup['mm volume']):
        assert np.isnan(vol) or vol >= 0.0, msg.format(i)
    for i,vol in enumerate(df_setup['water volume']):
        assert np.isnan(vol) or vol >= 0.0, msg.format(i)
    

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

def pip_mastermix(df_setup, outFH, src_labware):
    """Writing worklist commands for aliquoting mastermix.
    Using 1-asp-multi-disp with certain tips.
    """
    # split by mastermix names (different master mixes)
    MM_names = np.unique(df_setup['mm name'])

    outFH.write('C;Master mix\n')    
    for i,name in enumerate(MM_names):
        df = df_setup.loc[df_setup['mm name'] == name]
        # filtering out 'nan' values
        df = df.loc[pd.notnull(df_setup['mm volume'])]
        df.index = range(df.shape[0])
        
        # determing how many multi-disp
        max_vol = max(df['mm volume'])
        if max_vol > 180:
            n_disp = int(np.floor(900 / max_vol))  # using 1000 ul tips
        else:
            n_disp= int(np.floor(180 / max_vol))   # using 200 ul tips

        # making multi-disp object
        MD = Fluent.multi_disp()
        MD.SrcRackLabel = src_labware
        MD.SrcPosition = i+1
        MD.DestRackLabel = df.dest_labware
        MD.DestPositions = df.dest_location
        MD.Volume = df['mm volume']
        MD.NoOfMultiDisp = n_disp
        # writing
        outFH.write(MD.cmd() + '\n')


def pip_samples(df_setup, outFH):
    """Commands for aliquoting samples into distination plate
    """
    outFH.write('C;Samples\n')
    # filtering 'nan' from table
    x = pd.notnull(df_setup['sample labware'])
    y = pd.notnull(df_setup['sample location'])
    df = df_setup.loc[x & y]
    df.index = range(df.shape[0])
    
    # for each Sample, write out asp/dispense commands
    for i in range(df.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = df.ix[i,'sample labware']
        asp.Position = df.ix[i,'sample location']
        asp.Volume = round(df.ix[i,'sample volume'], 2)
        asp.LiquidClass = 'Water Contact Wet Single'
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df.ix[i,'dest_labware']
        disp.Position = df.ix[i,'dest_location']
        disp.Volume = round(df.ix[i,'sample volume'], 2)
        disp.LiquidClass = 'Water Contact Wet Single'
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

def pip_water(df_setup, outFH, src_labware):
    """Writing worklist commands for aliquoting water
    Using single asp-disp.
    """

    # filtering 'nan' from table
    x = pd.notnull(df_setup['water volume'])
    df = df_setup.loc[x]
    df.index = range(df.shape[0])
            
    for i in range(df.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = src_labware
        asp.Position = 1
        asp.Volume = round(df.ix[i,'water volume'], 2)
        asp.LiquidClass = 'Water Contact Wet Single'
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df.ix[i,'dest_labware']
        disp.Position = df.ix[i,'dest_location']
        disp.Volume = round(df.ix[i,'water volume'], 2)
        disp.LiquidClass = 'Water Contact Wet Single'
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

def add_error(x, error_perc):
    if x is None:
        return None
    return x * (1.0 + error_perc / 100.0)

def write_report_line(outFH, subject, volume, round_digits=1, error_perc=None):
    if volume is None:
        v = 'NA'
    else:
        if error_perc is not None:
            volume = add_error(volume, error_perc)
        v = round(volume, round_digits)
    outFH.write('{}:\t{}\n'.format(subject, v))

def write_report(df_setup, outFH):
    """Writing a report on the qPCR setup
    """
    # calculating total volumes
    n_rxn = df_setup.shape[0]
    ## total mastermix
    total_mm_volume = np.sum(df_setup['mm volume'])
    ## total water
    total_water_volume = np.sum(df_setup['water volume'])

    # report
    # number of samples
    outFH.write('# PCR REPORT\n')
    outFH.write('Number of total rxns:\t{}\n'.format(n_rxn))
    ## raw total volumes
    outFH.write('# Total reagent volumes (ul)\n')
    write_report_line(outFH, 'Master Mix', total_mm_volume)
    write_report_line(outFH, 'Water', total_water_volume)
    ## with pipetting error
    outFH.write('# Total reagent volumes + 10% more (ul)\n')
    write_report_line(outFH, 'Master Mix', total_mm_volume, error_perc=10)
    write_report_line(outFH, 'Water', total_water_volume, error_perc=10)
    
    # samples
    outFH.write('')

    
def main(args=None):
    # Input
    if args is None:
        args = parse_args()
    check_args(args)

    # Load input table
    df_setup = load_setup(args.setup, 
                          file_format=args.format)
    check_df_setup(df_setup)

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
         df_setup = reorder_384well(df_setup, 'dest_location')
    
    # # Writing out gwl file
    gwl_file = args.prefix + '.gwl'
    with open(gwl_file, 'w') as gwlFH:
        ## Master mix(es)
        pip_mastermix(df_setup, outFH=gwlFH, src_labware='1x24 tube runner[001]')
        ## Samples
        pip_samples(df_setup, outFH=gwlFH)
        ## Water
        pip_water(df_setup, outFH=gwlFH, src_labware='100ml_1_1')

    # Creating report file
    report_file = args.prefix + '_report.txt'
    with open(report_file, 'w') as repFH:
        write_report(df_setup, outFH=repFH)
        

    # Create windows-line breaks formatted versions
    gwl_file_win = Utils.to_win(gwl_file)
    report_file_win = Utils.to_win(report_file)

    # end
    return(gwl_file, gwl_file_win, report_file, report_file_win)


# main
if __name__ == '__main__':
    pass


