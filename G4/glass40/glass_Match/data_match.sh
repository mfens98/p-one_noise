#!/usr/bin/bash

START=$1
END=$2

dir="/data/p-one/mens/G4/glass40/"

source /cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/geant4.sh
echo "Arguments: $1 $2"

for (( i = START; i <= END; i++ )); do
	printf "$i/$END\n"
	"$dir"build/scint Seed "$i" | sed -ne '/^START$/{:a' -e 'n;p;ba' -e '}' | awk '/Graphics/ {exit} {print}' > "$dir"glass_Match/glass_ms/ms$i.txt
done

