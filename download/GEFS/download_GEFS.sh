#!/bin/bash
######################################################################
# Filename:    download_GEFS.sh
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to download GEFS data
# https://registry.opendata.aws/noaa-gefs-pds/ (data link)
#
######################################################################

### Inputs
INIT_DATE=$1
DATA_NAME=$2
ENS=$3

### Set up paths
PATH_TO_OUT="/expanse/nfs/cw3e/cwp186/data/downloaded/GEFS/${INIT_DATE}/"
PATH_TO_DATA="s3://noaa-gefs-pds/gefs.${INIT_DATE}/00/atmos/${DATA_NAME}p5/"
# PATH_TO_DATA="s3://noaa-gefs-pds/gefs.${INIT_DATE}/00/${DATA_NAME}/"

lead_array=(
000
# 006
# 012
# 018
# 024
# 030
# 036
# 042
# 048
# 054
# 060
# 066
# 072
# 078
# 084
# 090
# 096
# 102
# 108
# 114
# 120
# 126
# 132
# 138
# 144
# 150
# 156
# 162
# 168
)
## outer loop - loop through ensemble members
for i in ${!lead_array[*]}
do
    LEAD="${lead_array[$i]}"
    INPUT_FNAME="${PATH_TO_DATA}${ENS}.t00z.${DATA_NAME}.0p50.f${LEAD}"
    # INPUT_FNAME="${PATH_TO_DATA}${ENS}.t00z.${DATA_NAME}f${LEAD}"
    OUTPUT_FNAME="${PATH_TO_OUT}${ENS}.t00z.${DATA_NAME}.0p50.f${LEAD}"
    # echo ${OUTPUT_FNAME}
    /expanse/nfs/cw3e/cwp140/aws/local/bin/aws s3 cp --region us-east-1 ${INPUT_FNAME} ${OUTPUT_FNAME} --no-sign-request

done
