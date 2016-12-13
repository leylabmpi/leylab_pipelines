#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_leylab_pipelines
----------------------------------

Tests for `leylab_pipelines` module.
"""
# import
## batteries
import os
import sys
import unittest
## 3rd party
import pandas as pd
## package
from leylab_pipelines import Map2Robot


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
class Test_Map2Robot(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # args
    def test_parser_dest(self):
        # 96-well
        args = Map2Robot.parse_args(['--desttype', '96-well',
                                     '--deststart', '97',
                                     'mapfile'])
        with self.assertRaises(ValueError):
            Map2Robot.check_args(args)
        # 384-well
        args = Map2Robot.parse_args(['--desttype', '384-well',
                                     '--deststart', '389',
                                     'mapfile'])
        with self.assertRaises(ValueError):
            Map2Robot.check_args(args)

    def test_parser_tube(self):
        # mastermix <1
        args = Map2Robot.parse_args(['--mmtube', '0', 'mapfile'])
        with self.assertRaises(ValueError):
            Map2Robot.check_args(args)
        # mastermix >24
        args = Map2Robot.parse_args(['--mmtube', '25', 'mapfile'])
        with self.assertRaises(ValueError):
            Map2Robot.check_args(args)

    # import
    def test_load_map_txt(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        args = Map2Robot.parse_args([mapfile])
        df = Map2Robot.map2df(args)
        self.assertTrue(isinstance(df, pd.DataFrame))

    def test_load_map_xls(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.xls')
        args = Map2Robot.parse_args([mapfile])
        df = Map2Robot.map2df(args)
        self.assertTrue(isinstance(df, pd.DataFrame))

    # adding destination
    def test_load_map_txt(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        args = Map2Robot.parse_args([mapfile])
        df_map = Map2Robot.map2df(args)
        df_map = Map2Robot.add_dest(df_map)
        self.assertTrue(isinstance(df_map, pd.DataFrame))


    def test_load_map_deststart(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        args = Map2Robot.parse_args([mapfile])
        df_map = Map2Robot.map2df(args)
        df_map = Map2Robot.add_dest(df_map, dest_start=49)
        self.assertTrue(isinstance(df_map, pd.DataFrame))
        loc_start = df_map.loc[0,'dest_location']
        self.assertEqual(loc_start, 49)
