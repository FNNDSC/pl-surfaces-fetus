#!/bin/bash -e

# Author: Jennings Zhang <jenni_zh@protonmail.com>
# Date: 23 Aug 2019
# expansion to pial surface has yet to be implemented.
# -g flag should be conditional dependent on sulcification index.

bin_folder=/neuro/users/jennings.zhang/workspace/bin
mcubes="nohup perl $bin_folder/marching_cubes_fetus.pl"
ffit="nohup perl $bin_folder/fit_subplate.pl"
smthd="python $bin_folder/smoothness.py"
vertd="python $bin_folder/edgy.py"
chamfer_script="bash $bin_folder/chamfer.sh"
diemesh="python $bin_folder/diemesh.py"
distortion="python $bin_folder/distortion_angles.py"
minccalc="minccalc -quiet -clobber -byte -unsigned -expr "

usage="usage: $0 [-l|-r] [-g] labels.mnc"

function show_help () {
  cat << HELP
$usage

Creates surfaces for the subplate zone.

The painted segmentation volume input should have these labels:
    3=white matter
    4=intermediate zone

options:
  -h    show this help message and exit
  -l    -left side
  -r    -right side
  -g    age is between 30-32 Gestational weeks
  -c    clobber, redo everything and overwrite existing files
HELP
}

old_brain="-slow"
young_brain="-small"
side=""
clobber=""

while getopts ":hlrgc" opt; do
  case $opt in
  h   ) show_help && exit 0 ;;
  l   ) side='-left' ;;
  r   ) side='-right' ;;
  g   ) old_brain='-g' && young_brain='' ;;
  c   ) clobber='clobber' ;;
  \?  ) echo "Invalid option: -$OPTARG\nRun $0 -h for help."
    exit 1 ;;
  esac
done
shift $((OPTIND-1))
labels=$1

if [ "${labels:(-4)}" != ".mnc" ]; then
  echo "Input labels must be a MINC file."
  exit 1
fi

if [ -z "$side" ]; then
  current_folder="${PWD,,}"
  if [[ $current_folder == *"left"* ]]; then
    side="-left"
  elif [[ $current_folder == *"right"* ]]; then
    side="-right"
  else
    echo "Must specify side: -l or -r"
    exit 1
  fi
  echo "Warning: -l/-r not specified."
  echo "Assuming side is \"$side\" from cwd."
fi

if [ ! -f $labels ]; then
  echo "$labels file not found"
  exit 1
fi

# ==============================================================================

#tmpdir=$(mktemp -d -t subplate-$(date +%Hh%M,%S)-XXXXXXXXX)
mkdir -p cache
mkdir -p qc
# a prefix option for all output filenames would be nice but idc
layer3=wm_81920.obj
layer4=iz_81920.obj

layer3_mask=cache/wm_mask.mnc
layer4_mask=cache/iz_mask.mnc
chamfer3=cache/layer3_chamfer.mnc
dist_file_layer3=qc/wm_dist.txt
smth_file_layer3=qc/wm_smth.txt
area_file_layer3=qc/wm_area.txt
dist_file_layer4=qc/iz_dist.txt
smth_file_layer4=qc/iz_smth.txt
area_file_layer4=qc/iz_area.txt
cubes_log=cache/wm_cubes.log
fit_log=cache/iz_fit.log
thickness_tlink=subplate_thickness.txt
thickness_tnear=qc/subplate_thickness_tnear.txt
thickness_discrepancy=qc/tlink_minus_tnear_discrepancy.txt
highlight_volume=cache/highlight_missing_subplate.mnc
diemask=diemask.txt
vertexmask=not_subplate_mask.txt
mid_surface=cache/mid.obj
raw_angles=qc/raw_distortion_angles.txt

if [ -f "$layer3" -a -z "$clobber" ]; then
  echo "$layer3 found, skipping white matter extraction."
else
  echo "Using marching cubes to extract the white matter surface..."
  echo "Copy and paste the command below to monitor progress."
  echo "watch -n 0.1 tail --lines 40 $(realpath $cubes_log)"
  set -x
  $minccalc "A[0]>2.5" $labels $layer3_mask
  $mcubes $side  $young_brain \
           $layer3_mask $layer3 > $cubes_log 2>&1

  $chamfer_script -c 0.0 $layer3_mask $chamfer3
  volume_object_evaluate -linear $chamfer3 $layer3 $dist_file_layer3
  $smthd $layer3 $smth_file_layer3
  depth_potential -area_simple $layer3 $area_file_layer3
  { set +x; } 2> /dev/null

  echo "Creating $vertexmask..."
  $diemesh -keep $highlight_volume -boundary 2 $labels $layer3 $vertexmask

  if [ -f $layer4 ]; then
    echo "$layer4 is outdated."
    rm -v $layer4
  fi
fi

# if [ -f $diemask ]; then
#   echo "Found $diemask"
# else
#   unknown_filename=(*diemask*.txt*)
#   if [ -f $unknown_filename ]; then # find diemask.txt.backup
#     echo "Found manually created vertex mask $unknown_filename"
#     chmod --verbose 0444 $unknown_filename
#     cp -v $unknown_filename $diemask
#   else
#     echo "Creating $diemask..."
#     python /neuro/users/jennings.zhang/workspace/bin/diemesh.py
#           -keep cache/die.mnc $labels $layer3 diemask.txt
#   fi
# fi

if [ -f "$layer4" -a -z "$clobber" ]; then
  echo "$layer4 found, skipping surface_fit."
else
  echo "Fitting $layer3 inwards to the intermediate zone..."
  echo "Copy and paste the command below to monitor progress."
  echo "watch -n 0.1 tail --lines 40 $(realpath $fit_log)"
  set -x
  $minccalc "A[0]>3.5" $labels $layer4_mask
  $ffit $old_brain -evaluate $dist_file_layer4 \
           $layer4_mask $layer3 $layer4 > $fit_log 2>&1
  $smthd $layer4 $smth_file_layer4
  depth_potential -area_simple $layer4 $area_file_layer4
  { set +x; } 2> /dev/null
fi

# don't use -fwhm 20
# -tlink seems to be the implied default method
set -x
cortical_thickness -tlink $layer4 $layer3 $thickness_tlink
cortical_thickness -tnear $layer3 $layer4 $thickness_tnear
vertstats_math -old_style_file -sub $thickness_tlink $thickness_tnear $thickness_discrepancy
$distortion -mid $mid_surface $layer4 $layer3 $raw_angles
{ set +x; } 2> /dev/null

cat << EOF
See the results:
Display -spectral $labels {wm,iz}_81920.obj
brain-view iz_81920.obj qc/iz_dist.txt
EOF
