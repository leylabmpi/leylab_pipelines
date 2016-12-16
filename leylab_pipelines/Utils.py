# -*- coding: utf-8 -*-

# import
from __future__ import print_function
import os
import sys
import numpy as np

# functions
def make_range(x, set_zero_index=False):
    """Making a range from a string of comma-delimited and hyphenated number strings.
    If x = 'all', None returned
    x : string
    set_zero_index : change from 1-index to 0-index
    Example: 1,2,3,4
    Example: 1,2,5-6
    """
    x = str(x)
    if x.lower() == 'all':
        return None
    y = str(x).split(',')

    z = []
    for i in y:
        j = [int(x) for x in i.split('-')]
        if len(j) > 1:
            z = z + list(range(j[0], j[1]+1))
        else:
            z = z + j

    if set_zero_index:
        z = [x-1 for x in z]
        if min(z) < 0:
            raise ValueError('Negative values from setting zero index')

    return z


def check_gwl(gwl_file):
    """Checking that gwl in correct format
    gwl_file : input file name or file handle
    """
    # file read
    try:
        gwl_file = open(gwl_file, 'r')
    except ValueError:
        pass
    # checks
    for line in gwl_file:
        line = line.rstrip()
        if line is '':
            msg = 'Empty lines not allowed in gwl files'
            raise ValueError(msg)
        if not line[0] in ('A', 'D', 'R', 'W', 'F', 'C', 'B'):
            msg = '"{}" not a valid command ID'
            raise ValueError(msg.format(line[0]))
    # file close
    gwl_file.close()



# main
if __name__ == '__main__':
    pass
