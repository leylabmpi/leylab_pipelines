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
    Convert a QIIME-formatted mapping file to a GWL file, which is used by the TECAN
    robot to conduct the NGS amplicon PCR prep (ie., combining MasterMix, primers, samples, etc).

    EXTRA COLUMNS in MAPPING FILE:
    * "TECAN_sample_labware" = The sample labware name on the robot worktable
    * "TECAN_sample_location" = The well or tube location (a number)
    * "TECAN_primer_labware" = The primer plate labware name on the robot worktable
    * "TECAN_primer_location" = The well location (1-96 or 1-384)
    * "TECAN_sample_rxn_volume" = The volume of sample to use per PCR (ul)

    CONTROLS:
    * For the positive & negative controls, include them in the mapping file.
    * For "sample_labware", use [TODO: what name?]

    NOTES:
    * All volumes are in ul
    * Plate well locations are 1 to n-wells; numbering by column
    * PicoGreen should be added to the MasterMix *prior* to loading on robot
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('mapfile', metavar='MapFile', type=str,
                         help='A QIIME-formatted mapping file with extra columns (see below)')
    groupIO.add_argument('--rows', type=str, default='all',
                         help='Which rows of the mapping file to use (eg., "all"=all rows; "1-48"=rows1-48; "1,3,5-6"=rows1+3+5+6)')
    groupIO.add_argument('--prefix', type=str, default='TECAN_NGS_amplicon',
                         help='Output file name prefix')

    ## destination plate
    dest = parser.add_argument_group('Destination plate')
    dest.add_argument('--desttype', type=str, default='96-well',
                      choices=['96-well','384-well'],
                      help='Destination plate labware type')
    dest.add_argument('--deststart', type=int, default=1,
                      help='Start well number on destination plate')
    dest.add_argument('--rxns', type=int, default=3,
                      help='Number of replicate PCRs per sample')
    dest.add_argument('--destlabware', type=str,
                      default='96-well:96 well[006],384-well:384 well[006]',
                      help='Choices for the destination labware name base on --desttype')

    ## MasterMix
    mm = parser.add_argument_group('Master Mix')
    mm.add_argument('--mmtube', type=int, default=1,
                        help='MasterMix tube number')
    mm.add_argument('--mmvolume', type=float, default=13.1,
                        help='MasterMix volume per PCR')

    ## Primers
    primers = parser.add_argument_group('primers')
    primers.add_argument('--fpvolume', type=float, default=1.0,
                        help='Forward primer volume per PCR')
    primers.add_argument('--rpvolume', type=float, default=1.0,
                        help='Reverse primer volume per PCR')
    primers.add_argument('--fptube', type=int, default=0,
                        help='Forward non-bacode primer tube number (0 = barcoded primer on a plate)')
    primers.add_argument('--rptube', type=int, default=0,
                        help='Reverse non-bacode primer tube number (0 = barcoded primer on a plate)')

    ## Misc
    misc = parser.add_argument_group('Misc')
    misc.add_argument('--pcrvolume', type=float, default=25.0,
                        help='Total volume per PCR')
    misc.add_argument('--errorperc', type=float, default=10.0,
                        help='Percent of extra total reagent volume to include')

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
    if args.desttype == '96-well':
        destlimit = 96
    elif args.desttype == '384-well':
        destlimit = 384
    if args.deststart < 1 or args.deststart > destlimit:
        msg = 'Destination start well # must be in range: 1-{}'
        raise ValueError(msg.format(destlimit))
    # destination labware
    args.destlabware = {x.split(':')[0]:x.split(':')[1] for x in args.destlabware.split(',')}
    # rows in mapping file
    args.rows = Utils.make_range(args.rows, set_zero_index=True)
    # tube number
    if args.mmtube < 1 or args.mmtube > 24:
        msg = '{} tube # must be in range: 1-24'
        raise ValueError(msg.format('MasterMix'))
    if args.fptube < 0 or args.fptube > 24:
        msg = '{} tube # must be in range: 1-24 (or 0 if no tube)'
        raise ValueError(msg.format('Forward primer'))
    if args.rptube < 0 or args.rptube > 24:
        msg = '{} tube # must be in range: 1-24 (or 0 if no tube)'
        raise ValueError(msg.format('Reverse primer'))
    # volumes
    ## mastermix
    if args.mmvolume > args.pcrvolume:
        msg = 'MasterMix volume > total PCR volume'
        raise ValueError(msg)
    if args.mmvolume / 2.0 > args.pcrvolume:
        msg = 'WARNING: MasterMix volume > half of PCR volume'
        print(msg, file=sys.stderr)
    ## primers?

def map2df(mapfile, row_select=None):
    """Loading a mapping file as a pandas dataframe
    """
    # load via pandas IO
    if mapfile.endswith('.txt') or mapfile.endswith('.csv'):
        df = pd.read_csv(mapfile, sep='\t')
    elif mapfile.endswith('.xls') or mapfile.endswith('.xlsx'):
        xls = pd.ExcelFile(mapfile)
        df = pd.read_excel(xls)
    else:
        raise ValueError('Mapping file not in usable format')

    # selecting particular rows
    if row_select is not None:
        df = df.iloc[row_select]

    # return
    return df

def check_df_map(df_map, args):
    """Assertions of df_map object formatting
    * Assumes `sample` field = 1st column
    """
    # checking columns
    req_cols = ['TECAN_sample_labware', 'TECAN_sample_location',
                'TECAN_primer_labware', 'TECAN_primer_location',
                'TECAN_sample_rxn_volume']

    msg = 'Required column "{}" not found'
    for req_col in req_cols:
        if req_col not in df_map.columns.values:
            raise ValueError(msg.format(req_col))

    # checking for unique samples
    if any(df_map.duplicated(df_map.columns.values[0])):
        msg = 'WARNING: duplicated sample values in the mapping file'
        print(msg, file=sys.stderr)
    # checking for unique barcodes
    if any(df_map.duplicated(df_map.columns.values[1])):
        msg = 'WARNING: duplicated barcodes in the mapping file'
        print(msg, file=sys.stderr)

    # checking sample volumes
    msg = 'WARNING: sample volume > mastermix volume'
    for sv in df_map['TECAN_sample_rxn_volume']:
        if sv > args.mmvolume:
            print(msg, file=sys.stderr)
        if sv < 0:
            raise ValueError('Sample volume < 0')

def add_dest(df_map, dest_labware_index, dest_type='96-well',
             dest_start=1, rxn_reps=3):
    """Setting destination locations for samples & primers.
    Making a new dataframe with:
      [sample, sample_rep, dest_labware, dest_location]
    * For each sample (i):
      * For each replicate (ii):
        * plate = destination plate type
        * well = i * (ii+1) + (ii+1) + start_offset
    Joining to df_map
    """
    dest_start= int(dest_start)
    try:
        dest_labware = dest_labware_index[dest_type]
    except KeyError:
        raise KeyError('Destination labware type not recognized')

    # init destination df
    sample_col = df_map.columns[0]
    cols = [sample_col, 'TECAN_pcr_rxn_rep',
            'TECAN_dest_labware', 'TECAN_dest_location']
    ncol = len(cols)
    nrow = df_map.shape[0] * rxn_reps
    if dest_type == '96-well' and nrow > 97 - dest_start:
        nrow = 97 - dest_start
    elif dest_type == '384-well' and nrow > 385 - dest_start:
        nrow = 385 - dest_start
    df_dest = pd.DataFrame(np.nan, index=range(nrow), columns=cols)

    # filling destination df
    for i,(sample,rep) in enumerate(product(df_map.ix[:,0], range(rxn_reps))):
        # dest location
        dest_location = i + dest_start
        msg = 'WARNING: Not enough wells for the number of samples'
        msg = msg + '. Truncating to max samples that will fit on the plate'
        if dest_type == '96-well' and dest_location > 96:
            print(msg, file=sys.stderr)
            break
        elif dest_type == '384-well' and dest_location > 384:
            print(msg, file=sys.stderr)
            break
            # adding values DF
        df_dest.iloc[i] = [sample, rep+1, dest_labware, dest_location]

    # df join (map + destination)
    df_j = pd.merge(df_map, df_dest, on=sample_col, how='inner')
    assert df_j.shape[0] == df_dest.shape[0], 'map-dest DF join error'
    assert df_j.shape[0] > 0, 'DF has len=0 after adding destinations'

    # return
    return df_j


def reorder_384well(df, reorder_col):
    """Reorder values so that the odd, then the even locations are
    transferred. This is faster for a 384-well plate
    df: pandas.DataFrame
    reorder_col: column name to reorder
    """
    df['TECAN_sort_IS_EVEN'] = [x % 2 == 0 for x in df[reorder_col]]
    df.sort_values(by=['TECAN_sort_IS_EVEN', reorder_col], inplace=True)
    df = df.drop('TECAN_sort_IS_EVEN', 1)
    return df


def pip_mastermix(df_map, outFH, mmvolume=13.1, mmtube=1):
    """Writing worklist commands for aliquoting mastermix.
    Using 1-asp-multi-disp with 200 ul tips.
    Method:
    * calc max multi-dispense for 200 ul tips & mmvolume
    * for 1:n_dispense
      * determine how many disp left for channel (every 8th)
      * if n_disp_left < n_disp: n_disp = n_disp_left
      * calc total volume: n_disp * mmvolume
    """
    mmtube = int(mmtube)
    outFH.write('C;MasterMix\n')
    # currently: tips used just once
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = 'micro15[{0:0>3}]'.format(mmtube) #'1x24 Eppendorf Tube Runner no Tubes[001]'
        asp.Position = mmtube
        #asp.RackType = '1x24 Eppendorf Tube Runner no Tubes'
        asp.Volume = mmvolume
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_map.ix[i,'TECAN_dest_labware']
        disp.Position = df_map.ix[i,'TECAN_dest_location']
        disp.Volume = mmvolume
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

#    # just 1 command needed
#    rd = Fluent.reagent_distribution()
#    # source
#    ## tube rack ID
#    rd.SrcRackLabel = '1x24 Eppendorf Tube Runner no Tubes[001]'
#    rd.SrcRackType = '1x24 Eppendorf Tube Runner no Tubes'
#    ## tube start/end
#    rd.SrcPosStart = mmtube
#    rd.SrcPosEnd = mmtube
#
#    # destination
#    ## plate ID
#    dest_labware = set(df_map['TECAN_dest_labware'])
#    assert len(dest_labware) == 1, 'multiple destination labware not allowed!'
#    rd.DestRackLabel = list(dest_labware)[0]
#    ## destination start/end
#    rd.DestPosStart = int(min(df_map['TECAN_dest_location']))
#    rd.DestPosEnd = int(max(df_map['TECAN_dest_location']))
#
#    # other
#    rd.Volume = mmvolume
#
#    # write
#    outFH.write(rd.cmd() + '\n')
#    outFH.write('W;\n')


def pip_nonbarcode_primer(df_map, outFH, volume, tube):
    """Pipetting primers from tube.
    Assuming primer is aliquoted to all samples
    df_map : mapping file dataframe
    volume : volume to aliquot to each reaction
    tube : tube number (rackID assumed)
    """
    tube = int(tube)
    outFH.write('C;Non-barcoded primers\n')
    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = 'micro15[{0:0>3}]'.format(tube)
        asp.Position = tube
        asp.Volume = volume
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_map.ix[i,'TECAN_dest_labware']
        disp.Position = df_map.ix[i,'TECAN_dest_location']
        disp.Volume = volume
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

def pip_primers(df_map, outFH, fp_volume=0, rp_volume=0,
                fp_tube=0, rp_tube=0):
    """Commands for aliquoting primers
    """
    outFH.write('C;Primers\n')
    primer_plate_volume = 0
    # pipetting non-barcoded primers
    ## forward primer
    if fp_tube > 0 and fp_volume > 0:
        pip_nonbarcode_primer(df_map, outFH, fp_volume, fp_tube)
    else:
        primer_plate_volume += fp_volume
    ## reverse primer
    if rp_tube > 0 and rp_volume > 0:
        pip_nonbarcode_primer(df_map, outFH, rp_volume, rp_tube)
    else:
        primer_plate_volume += rp_volume

    # pipetting barcoded primers
    outFH.write('C;Barcoded primers\n')
    ## for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = df_map.ix[i,'TECAN_primer_labware']
        asp.Position = df_map.ix[i,'TECAN_primer_location']
        asp.Volume = primer_plate_volume
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_map.ix[i,'TECAN_dest_labware']
        disp.Position = df_map.ix[i,'TECAN_dest_location']
        disp.Volume = primer_plate_volume
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

def pip_samples(df_map, outFH):
    """Commands for aliquoting samples to each PCR rxn
    """
    outFH.write('C;Samples\n')
    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = df_map.ix[i,'TECAN_sample_labware']
        asp.Position = df_map.ix[i,'TECAN_sample_location']
        asp.Volume = df_map.ix[i,'TECAN_sample_rxn_volume']
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_map.ix[i,'TECAN_dest_labware']
        disp.Position = df_map.ix[i,'TECAN_dest_location']
        disp.Volume = df_map.ix[i,'TECAN_sample_rxn_volume']
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

def pip_water(df_map, outFH, pcr_volume=25.0,
              mm_volume=13.1, fp_volume=2.0, rp_volume=2.0):
    """Commands for aliquoting water to each PCR rxn
    """
    outFH.write('C;Water\n')
    # calculate the amount of water
    water_volume = []
    for i in range(df_map.shape[0]):
        samp_volume = df_map.ix[i,'TECAN_sample_rxn_volume']
        w_need = pcr_volume - (samp_volume + mm_volume + fp_volume + rp_volume)
        assert w_need >= 0, 'Water volume is negative: {}'.format(w_need)
        water_volume.append(w_need)
    df_map['TECAN_water_rxn_volume'] = water_volume

    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = df_map.ix[i,'TECAN_sample_labware']
        asp.Position = df_map.ix[i,'TECAN_sample_location']
        asp.Volume = water_volume[i]
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_map.ix[i,'TECAN_dest_labware']
        disp.Position = df_map.ix[i,'TECAN_dest_location']
        disp.Volume = water_volume[i]
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

def add_error(x, error_perc):
    if x is None:
        return None
    return x * (100 + error_perc)

def write_report_line(outFH, subject, volume, round_digits=2, error_perc=None):
    if volume is None:
        v = 'NA'
    else:
        if error_perc is not None:
            volume = add_error(volume, error_perc)
        v = round(volume, round_digits)
    outFH.write('{}:\t{}\n'.format(subject, v))

def write_report(df_map, outFH, pcr_volume, mm_volume,
                 fp_tube, fp_volume, rp_tube, rp_volume,
                 n_rxn_reps, error_perc=10.0):
    """Writing a report on
    """
    # calculating total volumes
    n_rxn = df_map.shape[0]
    ## total PCR
    total_pcr_volume = pcr_volume * n_rxn
    ## total mastermix
    total_mm_volume = pcr_volume * n_rxn
    ## total primer
    if fp_tube > 0 and fp_volume > 0:
        total_fp_volume = fp_volume * n_rxn
    else:
        total_fp_volume = None
    if rp_tube > 0 and rp_volume > 0:
        total_rp_volume = rp_volume * n_rxn
    else:
        total_rp_volume = None
    ## total water
    total_water_volume = sum(df_map['TECAN_water_rxn_volume'])

    # report
    # number of samples
    outFH.write('# PCR REPORT\n')
    outFH.write('Number of PCR replicates:\t{}\n'.format(n_rxn_reps))
    outFH.write('Number of total PCRs:\t{}\n'.format(n_rxn))
    ## rxn volumes
    outFH.write('# Volumes per PCR (ul)\n')
    write_report_line(outFH, 'Total', pcr_volume)
    write_report_line(outFH, 'Master Mix', mm_volume)
    write_report_line(outFH, 'Forward Primer', fp_volume)
    write_report_line(outFH, 'Reverse Primer', rp_volume)
    ## raw total volumes
    outFH.write('# Total volumes (ul)\n')
    write_report_line(outFH, 'Master Mix', total_mm_volume)
    write_report_line(outFH, 'Forward Primer', total_fp_volume)
    write_report_line(outFH, 'Reverse Primer', total_rp_volume)
    write_report_line(outFH, 'Water', total_water_volume)
    ## total volumes with error
    outFH.write('# Total volumes with (ul; {}% error)\n'.format(error_perc))
    write_report_line(outFH, 'Master Mix', total_mm_volume, error_perc=error_perc)
    write_report_line(outFH, 'Forward Primer', total_fp_volume, error_perc=error_perc)
    write_report_line(outFH, 'Reverse Primer', total_rp_volume, error_perc=error_perc)
    write_report_line(outFH, 'Water', total_water_volume, error_perc=error_perc)
    # samples
    outFH.write('')



def main(args=None):
    # Input
    if args is None:
        args = parse_args()
    check_args(args)
    # Import
    df_map = map2df(args.mapfile, row_select=args.rows)
    check_df_map(df_map, args)
    # Making destination dataframe
    df_map = add_dest(df_map, args.destlabware,
                      dest_type=args.desttype,
                      dest_start=args.deststart,
                      rxn_reps=args.rxns)
    # Reordering dest if plate type is 384-well
    if args.desttype == '384-well':
        df_map = reorder_384well(df_map, 'TECAN_dest_location')
    # GWL file construction
    ## gwl open
    gwl_file = args.prefix + '.gwl'
    gwlFH = open(gwl_file, 'w')
    ## mastermix
    pip_mastermix(df_map, gwlFH,
                  mmtube = args.mmtube,
                  mmvolume = args.mmvolume)
    ## primers
    pip_primers(df_map, gwlFH,
                fp_volume=args.fpvolume,
                rp_volume=args.rpvolume,
                fp_tube=args.fptube,
                rp_tube=args.rptube)
    ## samples
    pip_samples(df_map, gwlFH)
    ## water
    pip_water(df_map, gwlFH,
              pcr_volume=args.pcrvolume,
              mm_volume=args.mmvolume,
              fp_volume=args.fpvolume,
              rp_volume=args.rpvolume)
    ## gwl close
    gwlFH.close()
    # Report (total volumes; sample truncation; samples)
    report_file = args.prefix + '.report'
    with open(report_file, 'w') as repFH:
        write_report(df_map, outFH=repFH,
                     pcr_volume=args.pcrvolume,
                     mm_volume=args.mmvolume,
                     fp_tube=args.fptube,
                     fp_volume=args.fpvolume,
                     rp_tube=args.rptube,
                     rp_volume=args.rpvolume,
                     n_rxn_reps=args.rxns,
                     error_perc=args.errorperc)

    # Mapping file with destinations
    df_file = args.prefix + '_map.txt'
    df_map['TECAN_water_rxn_volume'] = df_map['TECAN_water_rxn_volume'].round(2)
    df_map['TECAN_dest_location'] = df_map['TECAN_dest_location'].astype(int)
    df_map['TECAN_pcr_rxn_rep'] = df_map['TECAN_pcr_rxn_rep'].astype(int)
    df_map.to_csv(df_file, sep='\t', index=False, na_rep='NA')

    # Return
    return (gwl_file, report_file, df_file)

# main
if __name__ == '__main__':
    pass




