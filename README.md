
# Pyshifts PyMOL plugin
PyShifts is a graphical analysis tool that utilize chemical shifts to assess the global quality of NMR structures of RNA. For more information on theoretical basis as well as our promising test results, see http://linktoManuscript . 
 
## Installation
Pyshifts is a plugin in PyMOL, an open source Python-enhanced molecular graphics tool. Python of version 2.7.10 and PyMOL is REQUIRED for Pyshifts.
#### 1. Python
- Python version of 2.7.10 (and 2.7.10 only), which is freely available at https://www.python.org/downloads/release/python-2710/.
- If your current Python version is not 2.7.10, and you prefer not to change it, you can use the following commands to create a Python 2.7.10 environment temporarily for Pyshifts:
        conda create -n pyshifts python=2.7.10
        source active pyshifts
After each use of Pyshifts, use 
        source deactivate pyshifts
to go back to your normal Python settings.

#### 2. PyMOL 
You can obtain PYMOL at sourceforge https://sourceforge.net/projects/pymol/.  

#### 3. Adding Pyshifts to PyMOL
- Download or clone this git repository.
- Open PyMOL and then go to Plugin -> Plugin manager -> Install new plugin, and choose the Pyshifts.py file in your local Pyshifts repository. For this step PyMOL need to be run with the Tcl/Tk interface, read more on PyMOL wiki https://pymolwiki.org/index.php/Plugins.

## Using Pyshifts
For more detailed instructions, watch video on youtube http://youtube.com/TOBEDONE
1. Load the object to be analyzed in PyMOL.

2. Run Pyshifts through Plugin -> Pyshifts. There are a progress bar, four tabs and one 'Exit' button in the pop-up window. The first tab 'Options' include basic options for Pyshifts and main 'Execute' button. In the second tab 'Error Analysis', table results will be shown and there are options about sorting tables. In the 'Advanced Options' tab, there are a few options on error offset and rendering error in PyMOL.

3. Following the simple steps below to analyze and visualize error
  - Change the 'PyMOL selection/ object' entry to the name of your target object and Execute. Pyshifts will automatically use both Larmord and Ramsey to predict chemical shifts of each states in your target object.
  - Change the 'Chemical Shift File' entry to the file path of experimental chemical shifts data of your target object.
  - Go to the second tab 'Error Analysis' and click on 'Compare shifts'. Error of each nuclei in each states will be shown in Error Table, and can also be visualized in PyMOL.

## Copyright Notice

The PyMOL Plugin source code in this file is copyrighted, but you can
freely use and copy it as long as you don't change or remove any of
the copyright notices.

----------------------------------------------------------------------
                      This PyMOL Plugin is Copyright (C) 2016 by 
           Jingru Xie <jingrux at umich dot edu> and Aaron T. Frank <afrankz at umich dot edu>

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
