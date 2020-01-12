
# install python packages in pymol python environment

from pymol import cmd

import pymol.externing

pymol.externing.conda("install pandas")
pymol.externing.conda("install scipy")
pymol.externing.conda("install -c anaconda scikit-learn")
pymol.externing.conda("install -c schrodinger pymol-psico")


