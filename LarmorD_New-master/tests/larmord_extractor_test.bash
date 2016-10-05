#!/bin/bash

ntest=1
cd tests
test=1
../bin/larmord_extractor -csfile measured_shifts_A038.dat -beta -6 -cutoff 10 A038.pdb > extractor_test.log
diff extractor_test.log extractor_test.dat > diff_extractor_test.log
lines=`wc -l diff_extractor_test.log | awk '{print $1}'`
if [[ $lines == 0 ]] 
then
  echo "TEST PASSED"
  rm extractor_test.log diff_extractor_test.log
else
  echo "TEST FAILED (see tests/extractor_test.log)"
fi

