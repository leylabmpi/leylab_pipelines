#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATADIR=$DIR/../../../tests/data

# small table datasets
echo "# Input: just accessions"
DB acc2taxID --tax $DATADIR/nucl_gss.acc2taxid.gz --no-header $DATADIR/accessions.txt
echo "# Input: table with accessions in 1 column"
DB acc2taxID $DATADIR/accessions_col3.txt --tax $DATADIR/nucl_gss.acc2taxid.gz 
echo "# multiple taxonomy db files"
DB acc2taxID $DATADIR/accessions_col3.txt --tax $DATADIR/nucl_gss.acc2taxid.gz $DATADIR/nucl_wg.acc2taxid


# default 
echo "# Input: just accessions"
#DB acc2taxID --type gss --no-header $DIR/../../../tests/data/accessions.txt


echo "# Input: table with accessions in 1 column"
#DB acc2taxID $DIR/../../../tests/data/accessions_col3.txt
