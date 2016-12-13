# -*- coding: utf-8 -*-

# import
from __future__ import print_function
import os
import sys
import pandas

from leylab_pipelines import Utils

class TECAN_Series(pandas.Series):
    @property
    def _constructor(self):
        return TECAN_Series

    def func_test(self):
        return 'Custom function'


class TECAN_DataFrame(pandas.DataFrame):
    """Pandas subclassed DataFrame for working with TECAN Fluent software
    """
    def __init__(self, *args, **kw):
        super(TECAN_DataFrame, self).__init__(*args, **kw)

    @property
    def _constructor(self):
        return TECAN_DataFrame

    _constructor_sliced = TECAN_Series

    def func_test(self):
        return 'Custom function'



# main
if __name__ == '__main__':
    pass
