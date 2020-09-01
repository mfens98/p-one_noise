#!/usr/bin/bash

ST=$1
EN=$2

(( ST = 1*ST ))
(( EN = 1*EN ))
echo "Arguments: $ST $EN"

for (( i = ST; i <= EN; i++ )); do
	printf "$i/$EN\n"
	/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3 /data/p-one/mens/G4/absorption/reanalyze/reanalyze.py /data/p-one/mens/G4/absorption/reanalyze/cedarLinks/"$i".txt
done
