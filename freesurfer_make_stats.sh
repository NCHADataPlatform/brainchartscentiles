#!/bin/bash

if [ -z "$1" ]
then
	echo "Make csv stats files from the Freesurfer aseg and cortical parcellation stats files"
  echo
  echo "$0 <subjects list>"
	exit
fi

# make aseg stats 
asegstats2table --subjectsfile=$1 --meas=volume --delimiter=comma --tablefile=stats_aseg.csv --skip --common-segs

# make cortical parcellation stats for aparc (DK), aparc.a200s (Destrieux), DKT (DKT) atlases
for h in lh rh ; do
  for m in area thickness volume ; do
   aparcstats2table --subjectsfile=$1 --meas=${m} --hemi=${h} --delimiter=comma --parc=aparc.DKTatlas --tablefile=stats_aparc.DKTatlas_${h}_${m}.csv --skip
   aparcstats2table --subjectsfile=$1 --meas=${m} --hemi=${h} --delimiter=comma --parc=aparc.a2009s --tablefile=stats_aparc.a2009s_${h}_${m}.csv --skip
   aparcstats2table --subjectsfile=$1 --meas=${m} --hemi=${h} --delimiter=comma --parc=aparc        --tablefile=stats_aparc_${h}_${m}.csv --skip
  done 
done
