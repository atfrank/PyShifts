#!/bin/bash

# extract features for training and testing sets
rm -f  train_features.txt
i=1
rnas="1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U "
for rna in $rnas
do
	[[ $i == 1 ]] && bash extract_features.bash $rna  | tee train_features.txt 
	[[ $i != 1 ]] && bash extract_features.bash $rna | grep -v frame | tee -a train_features.txt	
	i=999
done

i=1
rm -f  test_features.txt
rnas="1Z2J 2JWV 2K63 2K64 2K65 2K66 2LI4 2LPS 2LQZ 2LUN 2LX1 2M12 2M21 2M22 2M23 2M24 2M57 2M8K 2MEQ 2MHI 2MI0 2MIS 1SCL 1XHP 1ZC5 2QH2 2QH3 2QH4"
for rna in $rnas
do
	[[ $i == 1 ]] && bash extract_features.bash $rna | tee test_features.txt 
	[[ $i != 1 ]] && bash extract_features.bash $rna  | grep -v frame | tee -a test_features.txt
	i=999
done
