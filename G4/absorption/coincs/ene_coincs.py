#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import sys
sys.path.append('/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/python3.6/site-packages')
sys.path.append('/data/p-one/mens/modules')
import numpy as np
import coincidence as c

infile = sys.argv[1]

d = c.read_in_data(infile)

uptimes = c.get_index(d,2,'up')
uprad = c.get_index(d,4,'up')
upene = c.get_index(d,5,'up')

downtimes = c.get_index(d,2,'down')
downrad = c.get_index(d,4,'down')
downene = c.get_index(d,5,'down')

print("Up lengths: times: %d,  radii: %d, energies %d" %(len(uptimes),len(uprad),len(upene)))
print("Down lengths: times: %d,  radii: %d, energies %d" %(len(downtimes),len(downrad),len(downene)))


#find all coincidences first
print("Finding coincidencs")
all_coincs = c.find_sim_coincidence(uptimes,downtimes,uprad,downrad,upene,downene)

a = np.shape(all_coincs)
print(a)
try:
	if a[1]==6:
		print("Writing out!")
		ival = str(infile).split('p')[-1].split('.')[0]
		pid = (ival-(ival%8))/8
		c.write_data(all_coincs,'/data/p-one/mens/G4/absorption/coincs/all_coincs/','all'+str(pid)+'.txt',append=True)
	else:
		print("WRONG SHAPE CHANGE IT IN SCRIPT!")
except IndexError:
	pass

