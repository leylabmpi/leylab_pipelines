# -*- coding: utf-8 -*-
# import
## batteries
from __future__ import print_function
import os
import sys
import argparse
## package
from leylab_pipelines import Map2Robot
from leylab_pipelines import Dilute
from leylab_pipelines import QPCR



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
