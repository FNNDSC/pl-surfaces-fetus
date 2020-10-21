#!/bin/bash
# Author: Jennings Zhang <jenni_zh@protonmail.com>
# Date: 23 Aug 2019

usage="usage: $0 input_labels.mnc patched_labels.mnc"

function show_help () {
  cat << HELP
$usage

Paints a one-voxel boundary of value 3 inside the outer surface
of the intermediate zone where the subplate is discontinuous,
connecting the subplate over the hippocampus and fixing errors
on the lateral surface while artificially connecting it over
the diencephalonic region.

options:
  -h           show this help message and exit
  -k mark.mnc  give a filename to keep the intermediate file
               which marks changed voxels.
HELP
}

if [[ $1 == *"-h"* ]]; then
  show_help && exit 0
fi

tmpdir=$(mktemp -d -t patch-subplate-XXXX)
sp=$tmpdir/sp.mnc
sp_grow=$tmpdir/sp_grow.mnc
iz=$tmpdir/iz.mnc
iz_shrink=$tmpdir/iz_shrink.mnc
marks=$tmpdir/marks.mnc

while getopts ":hk:" opt; do
  case $opt in
  h   ) show_help && exit 0 ;;
  k   ) marks=$OPTARG ;;
  \?  ) echo "Invalid option: -$OPTARG\nRun $0 -h for help."
    exit 1 ;;
  esac
done
shift $((OPTIND-1))
# input=$1
# output=$2

if [ -z "$2" ]; then
  printf "%s\n%s\n" "$usage" "Missing filenames"
  exit 1
fi

if [ ! -f "$1" ]; then
  echo "$1 does not exist."
  exit 1
fi

set -e # quit on error

# dilate sp region inwards into the iz
minccalc -quiet -unsigned -byte -expr "A[0]>2.5&&A[0]<3.5" $1 $sp
# magic command that makes dilate_volume work
# dilate_volume from minc-toolkit-1.9.16-20180117-Ubuntu_16.04-x86_64
# is bugged and will label the dilated regions with 0.0039216
# it is not necessary if using
# CIVET/quarantines/Linux-x86_64/bin/dilate_volume (libminc version 2.1.02)
# mincreshape -quiet -image_range 0 255 $sp $rsp
dilate_volume $sp $sp_grow 1 6 1 > /dev/null

minccalc -quiet -unsigned -byte -expr "A[0]>3.5" $1 $iz
dilate_volume $iz $iz_shrink 0 6 1 > /dev/null

# to create a volume
#minccalc -quiet -unsigned -byte -expr
#"if(A[2]>3.5&&A[2]<4.5){if(A[1]>0.5){4}else if(A[0]<0.5){3}else{4}}else{A[2]}"
#$iz_shrink $sp_grow $labels $output
minccalc -quiet -unsigned -byte -expr "A[2]>3.5&&A[2]<4.5&&A[1]<0.5&&A[0]<0.5" \
          $iz_shrink $sp_grow $1 $marks

rm -r $tmpdir

minccalc -quiet -clob -unsigned -byte -expr "if(A[1]>0.5){3}else{A[0]}" $1 $marks $2
