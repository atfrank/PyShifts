
# Pyshifts
Pyshifts: A Pymol Plugin to Chemical Shift-Based Analysis of Biomolecular Ensembles

## Installation
Pyshifts is a plugin in PyMOL. Tested on Pymol >= v2.0 and Python3. Should work with Pymol < v2.0.

#### 1. Set up pyshifts environment
        conda create -n pyshifts
        source activate pyshifts
        conda install pandas
        conda install scipy
        conda install -c anaconda scikit-learn 
        source deactivate pyshifts
        
#### 2. PyMOL 
You can obtain PYMOL [here](https://pymol.org/2/).

#### 3. Adding Pyshifts to PyMOL
- Download or clone this git repository.
- Open PyMOL and then go to Plugin -> Plugin manager -> Install new plugin, and choose the Pyshifts.py file in your local Pyshifts repository. For this step PyMOL need to be run with the Tcl/Tk interface, read more on [PyMOL wiki](https://pymolwiki.org/index.php/Plugins).

#### 4. Get Lamord package (for computing non-exchangeable RNA chemical shifts)
- LarmorD can be obtained [here](https://umich.flintbox.com/?embed=true#technologies/bae4aa7e-2fbd-47dd-8886-533c147d75df) and it is free of charge if not for commercial use. 
- You also have to set `LARMORD_BIN` path in your environment. For example:

        export LARMORD_BIN=/Software/LarmorD/bin
        export PATH="${LARMORD_BIN}:$PATH"

#### 4. Get LarmorCα package (for computing protein backbone chemical shifts)
- LarmorCα can be obtained [here](https://github.com/atfrank/LARMORCA) and it is free of charge if not for commercial use. 
- You also have to set `LARMORCA_BIN` path in your environment. For example:

        export LARMORCA_BIN=/Software/LARMORCA/bin
        export PATH="${LARMORCA_BIN}:$PATH"

#### 5. Get BME package
- Install [this](https://github.com/KULL-Centre/BME) Bayesian Maximum Entropy (BME) library.
- Remember to add the library to your PYTHONPATH. For example:

        export BME=/home/XXX/GitHub/BME/
        export PYTHONPATH="${BME}:$PYTHONPATH"


#### 6. Get Psico library
- [optional] Install the [Pymol ScrIpt COllection (PSICO)](https://github.com/speleo3/pymol-psico). Improves performance of PyShifts when computing chemical shifts using LARMORD.
- Remember to add the library to your PYTHONPATH. For example:

        export PSICO=/home/XXX/GitHub/pymol-psico/
        export PYTHONPATH="${PSICO}:$PYTHONPATH"

## Using Pyshifts

1. Load the object to be analyzed in PyMOL.

2. Run Pyshifts through Plugin -> Pyshifts. There are a progress bar, four tabs and one 'Exit' button in the pop-up window. The first tab 'Options' include basic options for Pyshifts and main 'Predict or Load Chemical Shifts' button. In the second tab 'Error Analysis', table results will be shown and there are options about sorting tables. In the 'Advanced Options' tab, there are a few options on error offset and rendering error in PyMOL.

3. Following the simple steps below to analyze and visualize error
  - Change the 'PyMOL selection/ object' entry to the name of your target object and 'Predict or Load Chemical Shifts'. Pyshifts will automatically use both Larmord to predict chemical shifts of each states in your target object.
  - Go to the second tab 'Error Analysis' and click on 'Compare shifts'. Error of each nuclei in each states will be shown in Error Table, and can also be visualized in PyMOL.

## Copyright Notice

The PyMOL Plugin source code in this file is copyrighted, but you can
freely use and copy it as long as you don't change or remove any of
the copyright notices.
                      This PyMOL Plugin is Copyright (C) 2016 by 
           Jingru Xie <jingrux at umich dot edu>, Kexin Zhang <kexin at umich dot edu> and Aaron T. Frank <afrankz at umich dot edu>
                              All Rights Reserved

Permission to use, copy, modify, distribute, and distribute modified
versions of this software and its documentation for any purpose and
without fee is hereby granted, provided that the above copyright
notice appear in all copies and that both the copyright notice and
this permission notice appear in supporting documentation, and that
the name(s) of the author(s) not be used in advertising or publicity
pertaining to distribution of the software without specific, written
prior permission.

THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
