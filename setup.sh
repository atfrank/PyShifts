#!/bin/bash

# pre-requesites: conda, xQuartz

# create conda environment for pyshifts and set path
conda create -n pyshifts
conda activate pyshifts

#conda install -c schrodinger pymol=2.2
conda install  -y -c schrodinger pymol
conda install  -y -c anaconda -c schrodinger pandas
conda install  -y -c anaconda -c schrodinger pymol-psico
conda install  -y -c anaconda -c schrodinger scikit-learn
conda install  -y -c anaconda -c schrodinger scipy

echo "export PYSHIFTS_PATH=$(pwd)" >> ~/.bashrc
echo "export PATH=\$PYSHIFTS_PATH:\$PATH" >> ~/.bashrc

# get Larmord & LarmorCa dependencies
if ! [ "$(bash -c 'echo ${LARMORD_BIN}')" ]
then
    git clone --depth=1 https://github.com/karoka/LarmorD_New.git
    cd LarmorD_new
    make clean
    make
    echo "export LARMORD_BIN=$(pwd)/bin" >> ~/.bashrc
    echo "export PATH=\$LARMORD_BIN:\$PATH" >> ~/.bashrc
    cd ..
fi

if ! [ "$(bash -c 'echo ${LARMORCA_BIN}')" ]
then
    git clone --depth=1 https://github.com/atfrank/LARMORCA.git
    cd LARMORCA
    make clean
    make
    echo "export LARMORCA_BIN=$(pwd)/bin" >> ~/.bashrc
    echo "export PATH=\$LARMORCA_BIN:\$PATH" >> ~/.bashrc
    cd ..
fi

source ~/.bashrc
# Install python packages 
# If pymol installed from installer, then python packages have to be installed
# inside pymol
# If pymol installed from conda, then there is no conda inside pymol
# and packages can be directly installed with conda install

# pymol -cq pymol_setup.py > error_cath
# if ! grep -q "No module named 'conda'" error_cath
# then
#     echo "Python Packages Installed inside PyMOL."
# else
#     # conda env update --file conda_setup.yml 
#     echo "Python Packages Installed with conda."
# fi    
# rm error_cath


