#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import sys
sys.path.append('/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/python3.6/site-packages')
import pandas as pd 
import numpy as np 

#read in the datafile
def read_in_data(filepath):
    data = pd.read_csv(filepath,delimiter='\t',header=None)
    channel = (np.array(data[0]))
    edge = np.array(data[1])
    time = np.array(data[2])
    return np.array([channel,edge,time])

#parse the data for times where the noise rate is less than 20mHz
def parse(times):
	
	maxcount = 20 #hit per ms
	final_ms = np.ceil(max(times)/1e6) #find the final ms of the data so that our bins actually match the ms

	masterlist = []
	for i in range(int(final_ms)): # find how many bins are in a particular ms of a data file
		ms_bin = times[np.where(np.logical_and(times>i*1e6,times<(i+1)*1e6))]

		if len(ms_bin)<maxcount: #if less counts than our chosen maximum, keep these datapoints otherwise discard
			masterlist = masterlist + list(ms_bin)

	return np.array(masterlist) #return the array of kept points

def write_out(filepath,up,down): #write out the data to a file (file paths likely broken)
	fname = filepath.split('/')[-1]

	np.savetxt(filepath[:-(len(fname)+len("datafiles/"))]+"parsed/"+str(fname)+"_up_parsed.txt", up,delimiter="\t")
	np.savetxt(filepath[:-(len(fname)+len("datafiles/"))]+"parsed/"+str(fname)+"_down_parsed.txt", down,delimiter="\t")




filepath = sys.argv[1] #get the path of the datafile

data = read_in_data(filepath) #read in the data

#get the times for the up and down facing pmts (on the lowest threshold)
uptimes = data[2][np.where(np.logical_and(data[0]==1,data[1]==0))]
downtimes = data[2][np.where(np.logical_and(data[0]==5,data[1]==0))]

#run the parser over the up and down pmts
parsed_up = parse(uptimes)
parsed_down = parse(downtimes)

#write out the pared data arrays for both the up and down pmt
write_out(filepath,parsed_up,parsed_down)

