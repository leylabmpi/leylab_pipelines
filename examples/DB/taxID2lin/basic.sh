#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# default 
echo "# Running basic script and writing output to temporary directory"
DB taxID2lin -p 4 -o /tmp/NCBI_taxID2lin.txt

