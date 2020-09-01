#!/usr/bin/bash

hugearr4=(/data/p-one/gaertner/2007_heimdall/*DARK*)
filenum=0
for f in "${hugearr4[@]}"; do
    ln -s $f "newData/file.$filenum"
    let "filenum++"
done

