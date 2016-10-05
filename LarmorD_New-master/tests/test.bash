cd tests
../bin/larmord  -csfile measured_shifts_A003.dat \
-parmfile ../data/larmorD_alphas_betas.dat \
-reffile ../data/larmorD_reference_shifts.dat \
-accfile ../data/larmorD_accuracy_weight.dat \
-cutoff 15.0 A003.pdb
echo " "
../bin/larmord  -csfile measured_shifts_2ADT.dat \
-parmfile ../data/larmorD_alphas_betas_rna.dat \
-reffile ../data/larmorD_reference_shifts_rna.dat \
-accfile ../data/larmorD_accuracy_weight_rna.dat \
-cutoff 9999.0 2ADT.pdb
echo " "
../bin/larmord  -csfile measured_shifts_2ADT.dat \
-parmfile ../data/larmorD_alphas_betas_rna.dat \
-reffile ../data/larmorD_reference_shifts_rna.dat \
-accfile ../data/larmorD_accuracy_weight_rna.dat \
-printError \
-cutoff 9999.0 2ADT.pdb

echo "2LX1 test shifts"
../bin/larmord  -csfile measured_shifts_2LX1.dat \
-parmfile ../data/larmorD_alphas_betas_rna.dat \
-reffile ../data/larmorD_reference_shifts_rna.dat \
-accfile ../data/larmorD_accuracy_weight_rna.dat \
-cutoff 9999.0 2LX1.pdb

echo " "
echo "2LX1 test error"
../bin/larmord  -csfile measured_shifts_2LX1.dat \
-parmfile ../data/larmorD_alphas_betas_rna.dat \
-reffile ../data/larmorD_reference_shifts_rna.dat \
-accfile ../data/larmorD_accuracy_weight_rna.dat \
-printError \
-cutoff 9999.0 2LX1.pdb

echo " "
echo "A003 test error residue- and nucleus-based weights"
../bin/larmord  -csfile measured_shifts_A003.dat \
-parmfile ../data/larmorD_proteins_alphas_betas_cutoff_5.dat \
-reffile ../data/larmorD_proteins_reference_shifts_cutoff_5.dat \
-accfile ../data/larmorD_proteins_residue_accuracy_cutoff_5.dat \
-printError -residueBasedWeights \
-cutoff 5 A003.pdb

echo " "
echo "A003 test error nucleus-based weights"
../bin/larmord  -csfile measured_shifts_A003.dat \
-parmfile ../data/larmorD_proteins_alphas_betas_cutoff_5.dat \
-reffile ../data/larmorD_proteins_reference_shifts_cutoff_5.dat \
-accfile ../data/larmorD_proteins_accuracy_cutoff_5.dat \
-printError \
-cutoff 5 A003.pdb 

o " "
echo "A003 test error nucleus-based weights"
../bin/larmord  -csfile measured_shifts_A003.dat \
-parmfile ../data/larmorD_proteins_alphas_betas_cutoff_5.dat \
-reffile ../data/larmorD_proteins_reference_shifts_cutoff_5.dat \
-accfile ../data/larmorD_proteins_accuracy_cutoff_5.dat \
-corrfile correl_weights_zero.txt  \
-printError \
-cutoff 5 A003.pdb

