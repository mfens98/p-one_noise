#!/usr/bin/bash

cedarFiles=(datafiles/*.txt)
i=0
for f in "${cedarFiles[@]}"; do
	ln -s ../"$f" cedarLinks/"$i".txt
	(( i++ ))
done
