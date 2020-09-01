#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import sys
sys.path.append('/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/python3.6/site-packages')
import numpy as np
import pandas as pd
import os
import fcntl

def read_in_data(filepath):
    data = pd.read_csv(filepath,delimiter='\t',header=None)
    channel = (np.array(data[0]))
    edge = np.array(data[1])
    time = np.array(data[2])
    return np.array([channel,edge,time])

def parse(times):
    
    maxcount = 20 #ms
    final_ms = np.ceil(max(times)/1e6)

    masterlist = []
    for i in range(int(final_ms)):
        ms_bin = times[np.where(np.logical_and(times>i*1e6,times<(i+1)*1e6))]

        if len(ms_bin)<maxcount:
            masterlist = masterlist + list(ms_bin)

    return np.array(masterlist)

def find_coincidence(up, down):
    if len(up)>len(down):
        short = down;long = up; sdown = True
    else:
        short = up; long = down; sdown = False
    
    data_list = []
    print("Analyzing "+str(len(short))+" times...")
    for i in range(len(short)):
        t0 = short[i]
        tmin = t0-25;tmax = t0+25 #give 25ns for a 'coincidence'
        potential_list = long[np.where(np.logical_and(long>tmin,long<tmax))]
    
        if len(potential_list)>0: 
            for j in range(len(potential_list)):
                diff = t0-potential_list[j]
                if (sdown): data_list.append([diff,-diff]) #second element is always up-down
                else: data_list.append([diff,diff])
    
    data_list = np.array(data_list)
    return data_list

def write_data(data_arr,filepath):

    filename = filepath.split('/')[-1]
    
    np.savetxt("/data/p-one/mens/Coincidence_analysis/parsing/coincidences/stricter/processed/"+filename+"_coincs.txt",data_arr,delimiter="\t")

def write_sim_data(data_arr):

    with open("/data/p-one/mens/G4/real_coinc_run/coincidences_2.txt",'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        np.savetxt(f,data_arr,delimiter="\t")
        fcntl.flock(f, fcntl.LOCK_UN)



if sys.argv[1][:2]=='--':
    filepath = sys.argv[2]
    if sys.argv[1]=='--simulation':
        sim = True
    else: 
        raise NameError("improper argument: '"+str(sys.argv[1])+"'. Did you mean --simulation?")
        sys.exit(1)
else:
    filepath = sys.argv[1]
    sim=False


data = read_in_data(filepath)

uptimes = data[2][np.where(np.logical_and(data[0]==1,data[1]==0))]
downtimes = data[2][np.where(np.logical_and(data[0]==5,data[1]==0))]

#check type of data we are analyzing
if sim:
    coincs = find_coincidence(uptimes,downtimes)
    a =  np.shape(coincs)
    print("Shape of coincidence array is: ", a)
    try:
        if a[1]==2:
            print("writing data")
            write_sim_data(coincs)
    except:
        exit(0)
else:
    parsed_up = parse(uptimes)
    parsed_down = parse(downtimes)

    coincs = find_coincidence(parsed_up,parsed_down)

    write_data(coincs,filepath)


