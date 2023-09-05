#!/bin/bash

my_path="$(dirname -- ${BASH_SOURCE[0]})"
VisPos="/usr/bin/Rscript $my_path/VisPosterior.R"


if [ "$#" -lt 1 ]; then
    >&2 echo "Usage: $(basename $0) inFile.res [alpha = 0.05]";
    exit 1;
fi

inFile=$1
alpha=$2;
[ -z $alpha ] && alpha=0.05;

echo $inFile |
    grep -q '\.res$' || { >&2 echo "inFile must have the .res extension"; exit 1; }

baseDir=$(dirname $inFile);
baseName=$(basename $inFile ".res");
outName="$baseDir/$baseName.tiff";

$VisPos $inFile tiff $alpha

tiffFiles=($(find "$baseDir" -name "${baseName}_*.tiff" | sort))

#for f in "${tiffFiles[@]}"; do
#    echo $f;
#done
#echo "> $outName"

tiffcp "${tiffFiles[@]}" $outName

rm -f "${tiffFiles[@]}"
