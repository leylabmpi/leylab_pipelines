#!/usr/bin/env python

from __future__ import print_function

import os
import sys
from leylab_pipelines import Dilute



# main
if __name__ == '__main__':
    files = Dilute.main()
    for F in files:
        print('File written: {}'.format(F), file=sys.stderr)
