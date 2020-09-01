#!/usr/bin/bash

declare -a files
readarray -t files < /data/p-one/mens/Coincidence_analysis/parsing/pulse_check/checkFiles.txt

for f in "${files[@]}"; do
	echo "$f"
	/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3 /data/p-one/mens/Coincidence_analysis/parsing/pulse_check/get_coinc_times.py "$f"
done
