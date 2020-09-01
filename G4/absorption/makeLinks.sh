#!/bin/bash

ifiles=(coincs/datafiles/*.txt)

let i=0
for f in "${ifiles[@]}"; do
	ln -s /data/p-one/mens/G4/absorption/"$f" reanalyze/links/"$i".txt
	(( i++ ))
done

#cfiles=

#for f in "${cfiles[@]}"; do
#	ln -s "$f" reanalyze/links/"$i".txt
#done
