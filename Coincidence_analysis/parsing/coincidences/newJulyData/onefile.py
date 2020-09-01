#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import numpy as np
import sys
import os

files=sys.argv[1]

f_list = os.listdir('./'+str(files))

master = []
for file in f_list:
	data = np.loadtxt(files+file,delimiter="\t")
	if len(data.flatten()) == 2:
		data = list([data])
	if len(master)==0:
		master = list(data)
	else:
		master = master + list(data)

master = np.array(master)
print(np.shape(master))

np.savetxt('./'+str(files)+"_all.txt",master,delimiter="\t")
sys.exit(0)
