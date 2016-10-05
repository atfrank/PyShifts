# LarmorD: A Distance-Based Chemical Shift Predictor
 
- Currently predicts 1H, 13C and 15N RNA chemical shifts

## Install
```shell
$ cd /path/to/LarmorD/
$ make clean
$ make 
```

## Usage manual
```shell
$ bin/larmord -h
Usage: larmord [-options] <PDBfile>

Options:
         -csfile : chemical shifts list file (string)
         -parmfile : LarmorD parameter file (string)
         -trj : trajectory file (string)
         -skip : skip rate for reading trajectory (integer)
         -start : frame at which to start reading trajectory (integer)
         -stop : frame at which to stop reading trajectoryframe (integer)
         -identification : ID tag used in output (string)

```

## Examples
```shell
$ # predict chemical shifts from a coordinate file 
$ bin/larmord -csfile tests/measured_shifts.dat -identification 1SCL tests/file.pdb (PDB format)
$ # predict chemical shifts from a trajectory file (DCD format) 
$ bin/larmord -csfile tests/measured_shifts.dat -identification 1SCL -trj tests/file.dcd tests/file.pdb
```

## output
### format
_processor, trajectory-frame, residue-number, residue-name, nucleus, predicted-shifts, measured-shifts, random-coil-shifts, id-tag_

### example
```shell
$ bin/larmord -csfile tests/measured_shifts.dat -identification 1SCL tests/file.pdb
  
  0 1 1 GUA C4' 83.3042 80.854 83.1 1SCL
  0 1 1 GUA H4' 4.46852 4.533 4 1SCL
  0 1 1 GUA C1' 91.9503 89.246 93.1 1SCL
  0 1 1 GUA H1' 5.73867 5.826 5.05 1SCL
  0 1 1 GUA C8 139.384 136.761 137.9 1SCL
  0 1 1 GUA H8 7.88394 8.154 7.69 1SCL
  0 1 1 GUA C3' 74.5986 72.559 73.9 1SCL
  0 1 1 GUA H3' 4.58187 4.927 4.14 1SCL
  0 1 1 GUA C2' 75.7191 72.732 76.3 1SCL
  0 1 1 GUA H2' 4.60956 4.467 4.02 1SCL
  0 1 2 GUA C4' 82.4796 79.711 83.1 1SCL
  0 1 2 GUA H4' 4.4962 4.548 4 1SCL
  0 1 2 GUA C1' 92.4939 90.258 93.1 1SCL
  0 1 2 GUA H1' 5.7539 5.932 5.05 1SCL
  0 1 2 GUA C8 137.111 134.495 137.9 1SCL
  0 1 2 GUA H8 7.4922 7.556 7.69 1SCL
  0 1 2 GUA C3' 73.1824 70.475 73.9 1SCL
  0 1 2 GUA H3' 4.57466 4.559 4.14 1SCL
  0 1 2 GUA C2' 75.5983 72.954 76.3 1SCL
  0 1 2 GUA H2' 4.67613 4.718 4.02 1SCL
  ...
```

## Chemical shift list file
### file format
_residue-name, residue-number, nucleus, measured-shifts, error_

### example
```shell
$ head -10 tests/measured_shifts.dat
  
  CYT 6 C1' 91.438 0.00
  URA 7 C1' 92.388 0.00
  CYT 8 C1' 89.297 0.00
  ADE 9 C1' 85.504 0.00
  GUA 10 C1' 81.3 0.00
  GUA 1 C1' 89.246 0.00
  GUA 2 C1' 90.258 0.00
  GUA 3 C1' 90.657 0.00
  URA 4 C1' 90.835 0.00
  GUA 5 C1' 89.925 0.00
  ...
```

## LarmorD parameter file
### file format
_nucleus, neighbor-residue-name, neighbor-atom-name, alpha_

### example for H5''
```shell
$ head -68 data/parameters.txt

  H5'' GUA C1' 3.67815480447
  H5'' GUA C2' 1.63100734262
  H5'' GUA C3' 2.89281770448
  H5'' GUA C4' 2.3548429102
  H5'' GUA C5' 0.0321335101989
  H5'' GUA P -0.171617457218
  H5'' GUA O5' 1.57702037737
  H5'' GUA O3' 0.326168827858
  H5'' GUA C2 -6.64414148375
  H5'' GUA C4 -14.2514277865
  H5'' GUA C5 -18.4752206521
  H5'' GUA C6 -23.3660033205
  H5'' GUA C8 11.2256815971
  H5'' GUA N1 -15.0357687755
  H5'' GUA N2 18.8071249973
  H5'' GUA N3 7.49438758708
  H5'' GUA N7 6.47382380489
  H5'' GUA N9 -1.33876233081
  H5'' GUA O6 13.3592994749
  H5'' ADE C1' -6.87501637071
  H5'' ADE C2' -0.125715293229
  H5'' ADE C3' 0.934354550638
  H5'' ADE C4' 0.256912097145
  H5'' ADE C5' 0.571155032105
  H5'' ADE P -0.433798841887
  H5'' ADE O5' 1.5556382895
  H5'' ADE O3' 0.852536186953
  H5'' ADE C2 19.465835991
  H5'' ADE C4 -2.63457000853
  H5'' ADE C5 -6.20431462114
  H5'' ADE C6 -5.52564589385
  H5'' ADE C8 3.72779586352
  H5'' ADE N1 8.97415774656
  H5'' ADE N3 -3.65500648528
  H5'' ADE N6 -21.6733704635
  H5'' ADE N7 0.39668878159
  H5'' ADE N9 1.89834133117
  H5'' URA C1' -7.8361401694
  H5'' URA C2' -1.05475191362
  H5'' URA C3' -1.23351633568
  H5'' URA C4' 4.8561240664
  H5'' URA C5' 0.0
  H5'' URA P 1.86617611577
  H5'' URA O5' 1.01063675972
  H5'' URA O3' 0.949438715561
  H5'' URA C2 3.44987727435
  H5'' URA C4 -7.50237209988
  H5'' URA C5 -14.857005414
  H5'' URA C6 3.44656331225
  H5'' URA N1 2.71399930928
  H5'' URA N3 19.9073925365
  H5'' URA O4 23.4150612805
  H5'' CYT C1' 0.706640452452
  H5'' CYT C2' 0.166171606825
  H5'' CYT C3' -3.1035073294
  H5'' CYT C4' -1.57037444807
  H5'' CYT C5' 1.32626298398
  H5'' CYT P -0.440649244862
  H5'' CYT O5' -2.63138895134
  H5'' CYT O3' 1.92397549643
  H5'' CYT C2 8.20045630532
  H5'' CYT C4 1.29679400155
  H5'' CYT C5 4.77736668512
  H5'' CYT C6 1.05351000601
  H5'' CYT N1 2.78643038779
  H5'' CYT N3 6.38091065052
  H5'' CYT N4 -12.0144974986
  H5'' CYT O2 -2.71778957123
  ...
```

## License
```
  Copyright University of Michigan.
  This file is part of the LARMOR software suite and is made available under license.
  See: LICENSE_COMMERICAL and LICENSE_NON-COMMERICAL   for licensing info
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.  
  Author: Sean M. Law, Aaron T. Frank and Charles L. Brooks III

```
