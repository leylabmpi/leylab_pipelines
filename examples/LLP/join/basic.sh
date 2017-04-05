#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATADIR=$DIR/../../../tests/data

# default 
echo "# Input: small table join"
LLP join --join accession=accession --how left $DATADIR/accessions_col3.txt $DATADIR/nucl_wg.acc2taxid

echo "# Input: small table join with gzip file"
LLP join --join accession=accession --how left $DATADIR/accessions_col3.gz $DATADIR/nucl_wg.acc2taxid

echo "# Input: same table; inner join on 2 columns"
LLP join --join accession=accession,OTU=OTU --how left $DATADIR/accessions_col3.txt $DATADIR/accessions_col3.txt

echo "# Input: same table; inner join on 1 column with float dtype"
LLP join --join abundance=abundance --how left -L float -R float $DATADIR/accessions_col3.txt $DATADIR/accessions_col3.txt

echo "# Input: same table; inner join on 1 column with object dtype"
LLP join --join abundance=abundance --how left $DATADIR/accessions_col3.txt $DATADIR/accessions_col3.txt
