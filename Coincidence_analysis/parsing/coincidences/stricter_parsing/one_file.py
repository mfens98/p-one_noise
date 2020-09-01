#!/cvmfs/python3

import os
import numpy as np 

prefix = '/data/p-one/mens/Coincidence_analysis/parsing/coincidences/stricter/processed/'
files = os.listdir(prefix)
print(len(files))
masterlist = []
for file in files:
	data = np.loadtxt(prefix+file,delimiter="\t")
	if len(data.flatten()) == 2:
		data = list([data])
	if len(masterlist)==0:
		masterlist = list(data)
	else:
		masterlist = masterlist + list(data)

masterlist = np.array(masterlist)
np.savetxt('/data/p-one/mens/Coincidence_analysis/parsing/coincidences/stricter/master.gz', masterlist, delimiter="\t")
