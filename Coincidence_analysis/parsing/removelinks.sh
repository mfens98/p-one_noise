#!/usr/bin/bash

arr=(./datafiles/file.*)

for f in "${arr[@]}"; do

    unlink $f

done
