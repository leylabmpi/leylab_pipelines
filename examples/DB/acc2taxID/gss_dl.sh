#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATADIR=$DIR/../../../tests/data

# default 
echo "# Input: just accessions"
DB acc2taxID --type gss --no-header $DATADIR/accessions.txt


echo "# Input: table with accessions in 1 column"
#DB acc2taxID $DIR/../../../tests/data/accessions_col3.txt
