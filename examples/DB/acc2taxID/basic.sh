#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# default 
echo "# Input: just accessions"
#DB acc2taxID --type gss --no-header $DIR/../../../tests/data/accessions.txt

# pre-downloaded tax-dmp
DB acc2taxID --tax /tmp/nucl_gss.accession2taxid.gz \
    --no-header $DIR/../../../tests/data/accessions.txt

echo "# Input: table with accessions in 1 column"
#DB acc2taxID $DIR/../../../tests/data/accessions_col3.txt
