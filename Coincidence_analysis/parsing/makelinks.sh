#!/usr/bin/bash

hugearr=(/data/p-one/gaertner/2004_embla/unpacked/*DARK*)
hugearr2=(/data/p-one/gaertner/2005_fafnir/unpacked/*DARK*)
hugearr3=(/data/p-one/gaertner/2006_gungnir/*DARK*)
hugearr4=(/data/p-one/gaertner/2007_heimdall/*DARK*)
filenum=0
for f in "${hugearr[@]}"; do
    ln -s $f "datafiles/file.$filenum"

    let "filenum++"

done

for f in "${hugearr2[@]}"; do
    ln -s $f "datafiles/file.$filenum"
    let "filenum++"
done

for f in "${hugearr3[@]}"; do
    ln -s $f "datafiles/file.$filenum"
    let "filenum++"
done

for f in "${hugearr4[@]}"; do
    ln -s $f "datafiles/file.$filenum"
    let "filenum++"
done

