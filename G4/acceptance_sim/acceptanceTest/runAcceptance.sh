#!/bin/bash

rm -f "analysis_QE.txt"
printf "WARNING: This G4 simulation has not been modified since May 29, the setup may have changed\n"
END=180
for ((i=0;i<=END;i+=10)); do
    echo "$i"
    printf "$i\t" >> "analysis_QE.txt"
    ../K40accept n 1000000 angle "$i" Seed "$i" | grep -c up >> "analysis_QE.txt"

done

python3 makeplot.py
