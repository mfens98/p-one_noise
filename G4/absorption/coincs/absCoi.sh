#!/usr/bin/bash

START=$1
END=$2

dir="/data/p-one/mens/G4/absorption/"

source /cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/geant4.sh
echo "Arguments: $1 $2"

for (( i = START; i <= END; i++ )); do
	
	let "pid=$START/8"

	printf "$i/$END\n"
	"$dir"/build/K40coinc Seed "$i" > tmp2$i.txt
	parr=($(grep generated: tmp2$i.txt))
	n=${parr[2]}
	printf "$n\n" >> $dir/coincs/nGen/pid$pid.txt


	sed -ne '/^START$/{:a' -e 'n;p;ba' -e '}' tmp2$i.txt | awk '/Graphics/ {exit} {print}' > tmp$i.txt
	rm -f tmp2$i.txt
	
	cp tmp$i.txt "$dir"/coincs/datafiles/data$i.txt

	/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3 "$dir"/coincs/ene_coincs.py tmp$i.txt
	
	rm -f tmp$i.txt
done

