#!/bin/bash
MCI_A="test_3d.mci"
MCI_B="test_3d_B.mci"

rm -f 3d_cpu_RT*.mco 3d_cpu_Tonly*.mco *Ronly_A.mco *Ronly_B.mco

compile_and_run() {
    local MODE=$1
    local LABEL=$2
    
    echo "Compiling for ${LABEL} (PD_MODE=${MODE})..."
    gcc -O3 -DPD_MODE=${MODE} -o local_mcml_3d mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -lm
    
    echo "Running A format for ${LABEL}..."
    ./local_mcml_3d ${MCI_A} > /dev/null
    mv test_3d.mco 3d_cpu_${LABEL}_A.mco
    
    echo "Running B format for ${LABEL}..."
    ./local_mcml_3d ${MCI_B} > /dev/null
    mv test_3d_B.mco 3d_cpu_${LABEL}_B.mco
    
    echo "Done ${LABEL}."
}

compile_and_run 1 "Ronly"
compile_and_run 2 "Tonly"
compile_and_run 3 "RT"

echo ""
echo "==== Verification ==== "
ls -lh *.mco
