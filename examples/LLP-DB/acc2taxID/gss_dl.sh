#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATADIR=$DIR/../../../tests/data

# default 
echo "# Input: just accessions"
LLP-DB acc2taxID --type gss --no-header $DATADIR/accessions.txt

