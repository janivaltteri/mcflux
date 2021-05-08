#!/bin/bash

set -Eeuo pipefail

farray=()

for f in $1*.pdf
do
    farray+="$f "
done

pdfunite ${farray[*]} temp.pdf

mv temp.pdf $1-flux.pdf
