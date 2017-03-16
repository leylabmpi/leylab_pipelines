#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
qPCR_setup.py --prefix /tmp/TECAN_qPCR $DIR/../../tests/data/qPCR_setup/qPCR_Zach_plate1.xlsx

