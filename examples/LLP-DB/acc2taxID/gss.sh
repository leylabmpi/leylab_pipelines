#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATADIR=$DIR/../../../tests/data

# default 
echo "# Input: just accessions"
LLP-DB acc2taxID --tax /tmp/nucl_gss.accession2taxid.gz --no-header $DATADIR/accessions.txt

