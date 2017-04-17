#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATADIR=$DIR/../../../tests/data

# small table datasets
echo "# Input: just accessions"
LLP-DB convert $DATADIR/accessions.txt
