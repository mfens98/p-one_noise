#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import sys
sys.path.append('/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/python3.6/site-packages')
import pandas as pd
import os 
import numpy as np 

def read_in_data(filepath):
    data = pd.read_csv(filepath,delimiter='\t',header=None)
    channel = (np.array(data[0]))
    edge = np.array(data[1])
    time = np.array(data[2])
    return np.array([channel,edge,time])


def parse(times):
	
	maxcount = 20
 #ms
	final_ms = np.ceil(max(times)/1e6)

	count=0
	for i in range(int(final_ms)):
		ms_bin = times[np.logical_and(times>i*1e6,times<(i+1)*1e6)]

		if len(ms_bin)<maxcount and len(ms_bin)>0:
			count+=1

	return count

def write_out(values):

	print(min(values),"\t",max(values))


filepath = sys.argv[1]

data = read_in_data(filepath)

uptimes = data[2][np.where(np.logical_and(data[0]==1,data[1]==0))]
downtimes = data[2][np.where(np.logical_and(data[0]==5,data[1]==0))]

parsed_up = parse(uptimes)
parsed_down = parse(downtimes)

write_out([parsed_up,parsed_down])

