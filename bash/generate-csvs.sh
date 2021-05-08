#!/bin/bash

set -Eeuo pipefail

datadir=$1

mapfile -d $'\0' farray < <(find $1 -name "*.xlsx" -print0)

for i in "${farray[@]}"
do
    fstring=$i
    newext=".csv"
    fileout="${fstring/'.xlsx'/$newext}"
    echo $fileout
    xlsx2csv -s 2 -i --skipemptycolumns -t %H:%M:%S.0 $i $fileout
done

