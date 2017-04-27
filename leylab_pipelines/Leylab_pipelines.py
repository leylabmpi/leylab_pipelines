# -*- coding: utf-8 -*-
# import
## batteries
from __future__ import print_function
import os
import sys
import argparse
## package
### TECAN
from leylab_pipelines.TECAN import Map2Robot
from leylab_pipelines.TECAN import Dilute
from leylab_pipelines.TECAN import QPCR
### LLP-DB
from leylab_pipelines.DB import Convert
from leylab_pipelines.DB import Acc2TaxID
from leylab_pipelines.DB import TaxID2Lin
from leylab_pipelines.DB import TaxID2LinTbl
from leylab_pipelines.DB import EggNOG
### LLP
from leylab_pipelines import Join


def TECAN_arg_parse():
  desc = 'Tools for working with the TECAN robot' 
  # subcommand descriptions
  epi = 'SUBCOMMANDS:\n'
  epi = epi + '  map2robot - ' + Map2Robot.get_desc() + '\n'
  epi = epi + '  dilute - ' + Dilute.get_desc() + '\n'
  epi = epi + '  qPCR - ' + QPCR.get_desc() + '\n'
  
  # main command arg parser
  parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                   formatter_class=argparse.RawTextHelpFormatter)
  
  # subparsers
  subparsers = parser.add_subparsers(dest='subparser_name')
  parser_map2robot = Map2Robot.parse_args(subparsers=subparsers)
  parser_dilute = Dilute.parse_args(subparsers=subparsers)
  parser_QPCR = QPCR.parse_args(subparsers=subparsers)

  # parsing args
  args = parser.parse_args()
  
  # running function
  if args.subparser_name == 'map2robot':
    Map2Robot.main(args)
  elif args.subparser_name == 'dilute':
    Dilute.main(args)
  elif args.subparser_name == 'qPCR':
    QPCR.main(args)
  else:
    msg = 'Command not recognized: "{}"'
    raise ValueError(msg.format(args.subparser_name))


def LLP_DB_arg_parse():
  desc = 'Tools for working with public databases' 
  # subcommand descriptions
  epi = 'SUBCOMMANDS:\n'
  epi = epi + '  convert - ' + Convert.get_desc() + '\n'
  epi = epi + '  acc2taxID - ' + Acc2TaxID.get_desc() + '\n'
  epi = epi + '  taxID2lin - ' + TaxID2Lin.get_desc() + '\n'
  epi = epi + '  taxID2linTbl - ' + TaxID2LinTbl.get_desc() + '\n'
  epi = epi + '  eggnog - ' + EggNOG.get_desc() + '\n'
  
  # main command arg parser
  parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                   formatter_class=argparse.RawTextHelpFormatter)
  
  # subparsers
  subparsers = parser.add_subparsers(dest='subparser_name')
  parser_convert = Convert.parse_args(subparsers=subparsers)
  parser_taxID2lin = TaxID2Lin.parse_args(subparsers=subparsers)
  parser_taxID2lintbl = TaxID2LinTbl.parse_args(subparsers=subparsers)
  parser_acc2taxID = Acc2TaxID.parse_args(subparsers=subparsers)
  parser_eggnog = EggNOG.parse_args(subparsers=subparsers)
  # parsing args
  args = parser.parse_args()
  
  # running command
  if args.subparser_name.lower() == 'convert':
    Convert.main(args)
  elif args.subparser_name.lower() == 'taxid2lin':
    TaxID2Lin.main(args)
  elif args.subparser_name.lower() == 'taxid2lintbl':
    TaxID2LinTbl.main(args)
  elif args.subparser_name.lower() == 'acc2taxid':
    Acc2TaxID.main(args)
  elif args.subparser_name.lower() == 'eggnog':
    EggNOG.main(args)
  else:
    msg = 'Command not recognized: "{}"'
    raise ValueError(msg.format(args.subparser_name))


def LLP_arg_parse():
  desc = 'General bioinformatic tools'
  # subcommand descriptions
  epi = 'SUBCOMMANDS:\n'
  epi = epi + '  join - ' + Join.get_desc() + '\n'
  
  # main command arg parser
  parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                   formatter_class=argparse.RawTextHelpFormatter)
  
  # subparsers
  subparsers = parser.add_subparsers(dest='subparser_name')
  parser_join = Join.parse_args(subparsers=subparsers)
  # parsing args
  args = parser.parse_args()
  
  # running function
  if args.subparser_name == 'join':
    Join.main(args)
  else:
    msg = 'Command not recognized: "{}"'
    raise ValueError(msg.format(args.subparser_name))
