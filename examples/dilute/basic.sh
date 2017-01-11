#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
dilute.py --prefix /tmp/TECAN_dilute $DIR/../../tests/data/conc_file1.txt

