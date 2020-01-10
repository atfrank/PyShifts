
# Pyshifts
Pyshifts: A Pymol Plugin for Chemical Shift-Based Analysis of Biomolecular Ensembles

## Prerequisite
* [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

* [xQuartz](https://www.xquartz.org) (mac only)

* For first time PyMOL user, you will need a [PyMOL license file](https://pymol.org/2/buy.html?q=buy), as PyMOL is a commercial software.

## Quick Start

### Install Dependencies

```
git clone https://github.com/atfrank/PyShifts.git
cd Pyshifts
. setup.sh
```
The following command will invoke a pymol window
```
pymol
```


### Adding Pyshifts to PyMOL

In PyMOL window, go to `Plugin` -> `Plugin manager` -> `Install new plugin`, choose `Pyshifts.py` file in your local Pyshifts folder, and click `OK` on the next step. You will then see a pop-up message `Plugin "PyShiftsPlugin" has been installed`.


### Using Pyshifts


#### Load Object
Load the object to be analyzed in PyMOL, e.g. `2KOC_test.pdb` provided in `test/` folder, by typing `load test/2KOC_test.pdb` in pymol command line or dragging the file into PyMOL window.

#### Run Pyshifts

Run Pyshifts through `Plugin` -> `Legacy Plugins` ->   `Pyshifts`.


#### Analyze in Pyshifts
- Change the `PyMOL selection/ object` entry to the name of your target object, e.g. `2KOC_test` and clink on `Run`.

- Go to the second tab `Error Analysis`, and click on `Compare shifts`.

- Click on `Error table` or  `CS table` to save results.

## Advanced Options

### Perform BME Calculation
BME: a Bayesian Maximum Entropy (BME) library. For more detail please refer to the manuscript.

#### Installation
```
git clone https://github.com/KULL-Centre/BME
cd BME
echo "export BME=$(pwd)" >> ~/.bashrc
echo "export PYTHONPATH=\$BME:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

Once `BME` is installed, `Pyshifts` will be able to perform BME analysis and results will be shown in `Error table` ( `BME` column in `Error table` will be all `1.0` if `BME` not installed).

### Pyshifts Parameter Tuning

Pyshifts has four tabs and one `Exit` button. The first tab `Options` include basic options for `Pyshifts`, the second tab `Error Analysis` performs chemical shift comparison, displays table results and provides options in different ways of sorting.

The third tab `Advanced Options` contains functionality for parameter tuning. There are options on chemical shift error offset, accuracy, PyMOL rendering setting as well as machine learning clustering parameters.


## Publications
* `LarmorD`: Frank, Aaron T., Sean M. Law, and Charles L. Brooks III. "A simple and fast approach for predicting 1H and 13C chemical shifts: toward chemical shift-guided simulations of RNA." The Journal of Physical Chemistry B 118.42 (2014): 12168-12175.

* `LarmorCa`: Frank, Aaron T., et al. "Predicting Protein Backbone Chemical Shifts From Cα Coordinates: Extracting High Resolution Experimental Observables from Low Resolution Models." Journal of chemical theory and computation 11.1 (2014): 325-331.


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
