#!/usr/bin/bash

longarr=(/data/p-one/gaertner/2004_embla/unpacked/*DARK*)
count=0
for f in "${longarr[@]}"; do
    python3 make_pdfs.py "$f"
    let count+=1
    printf "$count/${#longarr[@]}\r"

done
