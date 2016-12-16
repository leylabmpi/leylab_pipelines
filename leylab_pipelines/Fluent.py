# -*- coding: utf-8 -*-

# import
from __future__ import print_function
import os
import sys
import numpy as np


def xstr(x):
    if x is None:
        return ''
    else:
        return x


class reagent_distribution():
    """Commands for aliquoting mastermix
    *AspirateParameters*
    SrcRackLabel
    SrcRackID
    SrcRackType
    SrcPosStart
    SrcPosEnd
    *DispenseParameters*
    DestRackLabel
    DestRackID
    DestRackType
    DestPosStart
    DestPosEnd
    *Other*
    Volume = How much volume to asp/disp?
    LiquidClass = Which liquid class to use (default: 'Waer Free Multi')?
    NoOfDiTiReuses = How many times to reuse tips?
    NoOfMultiDisp = How many multi-dispenses?
    Direction = Which way to pipette (default:0)?
    ExcludedDestWell = Semi-colon separated string of locations to exclude

    *WashParameters*
    None?

    # Example: R;100ml_2;;Trough 100ml;1;1;96 Well Skirted PCR[003];;96 Well Skirted PCR;1;96;20;Water Free Multi;1;5;0
    """

    def __init__(self, ):
        self._ID = 'R'
        # aspirate parameters
        self.SrcRackLabel = None
        self.SrcRackID = None
        self.SrcRackType = None
        self.SrcPosStart = 1
        self.SrcPosEnd = 1
        # dispense parameters
        self.DestRackLabel = None
        self.DestRackID = None
        self.DestRackType = None
        self.DestPosStart = 1
        self.DestPosEnd = 1
        # other
        self.Volume = 1
        self.LiquidClass = 'Water Free Multi'
        self.NoOfDiTiReuses = 1
        self.NoOfMultiDisp = 5
        self.Direction = 0
        self.ExcludedDestWell = None
        self.key_order = ['_ID',
                          'SrcRackLabel', 'SrcRackID', 'SrcRackType',
                          'SrcPosStart', 'SrcPosEnd',
                          'DestRackLabel', 'DestRackID', 'DestRackType',
                          'DestPosStart', 'DestPosEnd',
                          'Volume', 'LiquidClass', 'NoOfDiTiReuses',
                          'NoOfMultiDisp', 'Direction', 'ExcludedDestWell']

    def xstr(self, x):
        if x is None:
            return ''
        else:
            return x

    def cmd(self):
        # list of values in correct order
        vals = [getattr(self, x) for x in self.key_order]
        # None to blank string
        vals = [xstr(x) for x in vals]
        # convert all to strings
        vals = [str(x) for x in vals]
        # return
        return ';'.join(vals)


class asp_disp():
    """Commands for aliquoting mastermix
    *Parameters*
    RackLabel
    RackID
    RackType
    Position
    TubeID
    Volume
    LiquidClass
    TipType
    TipMask
    ForceRack
    MinDetected
    """

    def __init__(self):
        self._ID = ''
        # aspirate parameters
        self.RackLabel = None
        self.RackID = None
        self.RackType = None
        self._Position = 1
        self.TubeID = None
        self.Volume = None
        self.LiquidClass = 'Water Free Single'
        self.TipType = None
        self.TipMask = None
        self.key_order = ['_ID',
                          'RackLabel', 'RackID', 'RackType',
                          'Position', 'TubeID', 'Volume',
                          'LiquidClass', 'TipType', 'TipMask']

    def cmd(self):
        # list of values in correct order
        vals = [getattr(self, x) for x in self.key_order]
        # None to blank string
        vals = [xstr(x) for x in vals]
        # convert all to strings
        vals = [str(x) for x in vals]
        # return
        return ';'.join(vals)

    @property
    def Position(self):
        return self._Position

    @Position.setter
    def Position(self, value):
        self_Position = int(value)

class aspirate(asp_disp):
    def __init__(self):
        asp_disp.__init__(self)
        self._ID = 'A'


class dispense(asp_disp):
    def __init__(self):
        asp_disp.__init__(self)
        self._ID = 'D'



# main
if __name__ == '__main__':
    pass
