#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
map2robot.py --prefix /ebio/abt3_projects/databases/TECAN/worklists/basic_96well \
  $DIR/../../tests/data/basic_96well.txt

# 394-well source plate
map2robot.py --prefix /ebio/abt3_projects/databases/TECAN/worklists/basic_384well \
  $DIR/../../tests/data/basic_384well.txt

# 96-well source plate, 384-well destination plate
map2robot.py --prefix /ebio/abt3_projects/databases/TECAN/worklists/basic_96-384 \
  --desttype 384 \
  $DIR/../../tests/data/basic_96well.txt
