#!/bin/bash

date=`date`
echo "Starting AT "$date

lag=5
yyyy=`date -d '-'$lag' hours' -u +%Y`
mm=`date -d '-'$lag' hours' -u +%m`
dd=`date -d '-'$lag' hours' -u +%d`
hh=`date -d '-'$lag' hours' -u +%H`

PATH_TO_REPO=/data/projects/operations/ICT_IWT/

rm -f ${PATH_TO_REPO}figs/GFS_*_IntWest_latest_F*.png

cd ${PATH_TO_REPO}


filename="/data/projects/external_datasets/GFS/processed/"$yyyy$mm$dd$hh"/"$yyyy$mm$dd$hh"_F168.grb2"
while true; do
  if [[ -e "$filename" ]]; then                   # Check if the file exists
    filesize1=$(stat --format="%s" "$filename")   # Get the file size
    sleep 5                                      # Wait a few seconds
    filesize2=$(stat --format="%s" "$filename")   # Get the file size again

    if [[ "$filesize1" == "$filesize2" ]]; then  # Compare file sizes
      break                                      # Exit the loop if file size is not changing
    else
      echo $filename" is still being written (file size is changing)."
      sleep 30
    fi
  else
    echo $filename" does not exist."
    sleep 30
  fi
done

echo $filename" ready for processing"

date=`date`
echo "STARTING MAKING ICT Plots at "$date

apptainer exec -e --bind /data:/data /data/projects/operations/ICT_IWT/env/ICT_IWT.sif python /data/projects/operations/ICT_IWT/run_tool.py "GFS" "$yyyy$mm$dd$hh"

cd ${PATH_TO_FIGS}figs/
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih GFS_*_IntWest_latest_F*.png /data/projects/website/mirror/htdocs/Projects/CO_landfalling_ARs/images/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

date=`date`
echo "Finished at "$date


exit





