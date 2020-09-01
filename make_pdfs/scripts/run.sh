#!/bin/bash

./make_pdfs.sh
./make_pdfs2.sh
./make_pdfs3.sh

tar -czf png_files.tgz pdfs/ ylim_pdfs/

cp png_files.tgz /data/p-one/mens/make_pdfs/
