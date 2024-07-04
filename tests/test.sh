#!/bin/bash
BINDIR=../bin
pycmd=python3
TESTDIR=.

# echo "***TEST: LiH static test no frozen core real time propagation to check energy conservation"
# START_TIME=$SECONDS
# ${BINDIR}/mcend LiH.static.inp-R > lih-nofc-R.txt
# echo "closed-shell nofc time" $(($SECONDS - $START_TIME))"s"
# ${pycmd} check_energy_conserved.py expec.t
# 
# echo "***TEST: LiH static test no frozen core"
# START_TIME=$SECONDS
# ${BINDIR}/mcend LiH.static.inp > lih-nofc.txt
# echo "closed-shell nofc time" $(($SECONDS - $START_TIME))"s"
# mv expec.t expec.t-lih-cs-nofc
# echo "Comparing with the reference values"
# ${pycmd} compare_expec_v2.py expec.t-lih-cs-nofc expec.t-lih-cs-ref
# 
# START_TIME=$SECONDS
# echo "***TEST: LiH acf test no frozen core"
# ${BINDIR}/mcend LiH.acf.inp > lih-nofc.txt
# echo "closed-shell nofc time" $(($SECONDS - $START_TIME))"s"
# mv expec.t expec.t-lih-cs-acf-nofc
# echo "Comparing with the reference values"
# ${pycmd} compare_expec_v2.py expec.t-lih-cs-acf-nofc expec.t-lih-cs-acf-ref
# 
# echo "LiH static test closed shell frozen core"
# START_TIME=$SECONDS
# ${BINDIR}/mcend LiH.static.inp-cs > lih-cs.txt
# echo "closed-shell time" $(($SECONDS - $START_TIME))"s"
# mv expec.t expec.t-lih-cs
# echo "Comparing with the reference values"
# ${pycmd} compare_expec_v2.py expec.t-lih-cs expec.t-lih-cs-fc-ref
# 
# START_TIME=$SECONDS
# echo "LiH acf test closed shell frozen core"
# ${BINDIR}/mcend LiH.acf.inp-cs > lih-acf-cs.txt
# echo "acf time" $(($SECONDS - $START_TIME))"s"
# mv expec.t expec.t-lih-acf-cs
# echo "Comparing with the reference values"
# ${pycmd} compare_expec_v2.py expec.t-lih-acf-cs expec.t-lih-acf-cs-fc-ref
# 
# START_TIME=$SECONDS
# echo "LiH static test open shell"
# ${BINDIR}/mcend LiH.static.inp-os > lih-os.txt
# echo "open-shell time" $(($SECONDS - $START_TIME))"s"
# mv expec_spinorbital.t expec_spinorbital.t-lih-os
# echo "Comparing with the closed-shell values"
# ${pycmd} compare_expec_v2.py expec.t-lih-cs expec_spinorbital.t-lih-os
# 
# START_TIME=$SECONDS
# echo "LiH acf test open shell"
# ${BINDIR}/mcend LiH.acf.inp-os > lih-acf-os.txt
# mv expec_spinorbital.t expec_spinorbital.t-lih-acf-os
# echo "open-shell time" $(($SECONDS - $START_TIME))"s"
# ${pycmd} compare_expec_v2.py expec.t-lih-acf-cs expec_spinorbital.t-lih-acf-os
# 
START_TIME=$SECONDS
echo "LiH+ static test"
${BINDIR}/mcend LiH.static.inp-p > lih-p.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"
mv expec_spinorbital.t expec_spinorbital.t-lih-p
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec_spinorbital.t-lih-p expec_spinorbital.t-lih-p-ref
 
echo "H2 static test closed shell"
START_TIME=$SECONDS
${BINDIR}/mcend H2.static.inp-cs > h2-cs.txt
echo "closed-shell time" $(($SECONDS - $START_TIME))"s"
mv expec.t expec.t-h2-cs
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec.t-h2-cs expec.t-h2-cs-ref

START_TIME=$SECONDS
echo "H2 static test open shell"
${BINDIR}/mcend H2.static.inp-os > h2-os.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"
mv expec_spinorbital.t expec_spinorbital.t-h2-os
${pycmd} compare_expec_v2.py expec.t-h2-cs expec_spinorbital.t-h2-os

echo "OH static test"
START_TIME=$SECONDS
${BINDIR}/mcend OH.static.inp > oh.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"
mv expec_spinorbital.t expec_spinorbital.t-oh
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec_spinorbital.t-oh expec_spinorbital.t-oh-ref

# clean up after running the tests
rm *.txt *.chk *.t *.ij *.list *.t-h2-os *.t-lih-acf-os *.t-lih-os *.t-lih-p 
rm *.t-oh *.t-h2-cs *.t-lih-acf-cs *.t-lih-cs *.t-lih-cs-acf-nofc *.t-lih-cs-nofc 
rm finalpsi* finalrho *.Rt