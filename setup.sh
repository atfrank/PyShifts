# pre-requesites: conda, xQuartz
# create conda environment for pyshifts and set path
conda init bash
conda create -n pyshifts
conda activate pyshifts

# install Python dependenices
conda install -y -c schrodinger pymol=2.4
conda install -y -c anaconda -c schrodinger pandas
conda install -y -c anaconda -c schrodinger pymol-psico
conda install -y -c anaconda -c schrodinger scikit-learn
conda install -y -c anaconda wget

echo "# added during PYSHIFTS installation" >> ~/.bashrc
echo "export PYSHIFTS_PATH=$(pwd)" >> ~/.bashrc
echo "export PATH=\$PYSHIFTS_PATH:\$PATH" >> ~/.bashrc

# install BME
if ! [ "$(bash -c 'echo ${BME}')" ]
then    
    wget https://github.com/KULL-Centre/BME/archive/refs/tags/v1.0.tar.gz
    tar -xvf v1.0.tar.gz
    mv BME-1.0 BME
    cd BME
    echo "export BME=$(pwd)" >> ~/.bashrc
    echo "export PYTHONPATH=\$BME:\$PATH" >> ~/.bashrc
    cd ..
fi

# get Larmord and LarmorCa dependencies
if ! [ "$(bash -c 'echo ${LARMORD_BIN}')" ]
then
    git clone --depth=1 https://github.com/karoka/LarmorD_New.git
    cd LarmorD_New
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
