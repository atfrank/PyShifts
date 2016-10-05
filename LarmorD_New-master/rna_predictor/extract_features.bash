#!/bin/bash
if [[ $# -ne 1 ]]
then
	echo "  Description: extract LARMORD features"
	echo "  Usage: $0 < PDBIDs>"
	echo "  Usage: $0 2LX1"
else	
	PDB=$1
	
	# Setup Environment
	source ~/.bashrc
	larmorDHome=/Users/atfrank/Dropbox/Documents/GitSoftware/LarmorD_New/
	dataHome=/Users/atfrank/Dropbox/Documents/Assignment-Testing/Resolving_RNA_Structure_Using_Unassigned_Spectra/assignment_pools/simRNA/chemicalShiftData/
	
	# Get native ensemble
	mkdir -p ${PDB}	
	split_pdb.bash ${PDB} ${PDB}/model 0 &> ${PDB}/split.out
	models=`grep contains ${PDB}/split.out | awk '{print $3}'`	
	models=`seq 1 $models`	

  # Get data
	data=${dataHome}/measured_shifts_${PDB}.dat
	data_correct=${dataHome}/measured_shifts_corrected_${PDB}.dat  
	[[ -s ${data_correct} ]] && data=${data_correct}	
	cat $data | sort -n -k2 > ${PDB}/cs_1.txt
	data=${PDB}/cs_1.txt
  
  # Loop over NMR bundle and extract structure features
  for m in $models
  do
  	id=`echo "${PDB}_${m}"`
  	${larmorDHome}/bin/larmord_extractor -rna -ring -cutoffRing 10.0 -identification ${id} -csfile ${data} -beta -3 -cutoff 9999 ${PDB}/model_${m}.pdb  	
  done
fi
