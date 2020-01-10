#!/bin/bash

# create virtual environment venv
pip install virtualenv
virtualenv venv

# activate virtual env
source venv/bin/activate

# install python libraries
pip install pandas
pip install scipy
pip install -U scikit-learn

# download PyMOL from website

# set path for pymol and python
export PATH="/Applications/PyMOL.app/Contents/bin:$PATH" # PyMOL location
export PYTHONPATH=/Applications/PyMOL.app/Contents/lib/python3.7/site-packages/pmg_tk/startup:$PYTHONPATH
export PATH="/Users/kexin/Documents/local_software/pyshifts/venv/bin:$PATH" # virtual env location
export PYTHONPATH="/Users/kexin/Documents/local_software/pyshifts/venv/lib/python3.7/site-packages:$PYTHONPATH"

# copy Meter.py to PyShifts installation location
cp Meter.py /Applications/PyMOL.app/Contents/lib/python3.7/site-packages/pmg_tk/startup

# set path for LarmorD
export LARMORD_BIN=/Users/kexin/Documents/GitHub/LarmorD_New/bin
export PATH="${LARMORD_BIN}:$PATH"
# set path for LarmorCA
export LARMORCA_BIN=/Users/kexin/Documents/GitHub/LARMORCA/bin
export PATH="${LARMORCA_BIN}:$PATH"
# set path for BME
export BME=/Users/kexin/Documents/GitHub/BME/
export PYTHONPATH="${BME}:$PYTHONPATH"
# set path for psico
export PSICO=/Users/kexin/Documents/GitHub/pymol-psico/
export PYTHONPATH="${PSICO}:$PYTHONPATH"