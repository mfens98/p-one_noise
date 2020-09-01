#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import sys
sys.path.append('/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/python3.6/site-packages')
import numpy as np
import pandas as pd
import os

def read_in_data(filepath):
    data = pd.read_csv(filepath,delimiter='\t',header=None)
    channel = (np.array(data[0]))
    edge = np.array(data[1])
    time = np.array(data[2])
    return np.array([channel,edge,time])

def parse(times): #parse by the ms and if we have more than 20 hits throw out
    
    maxcount = 20
 #ms
    final_ms = np.ceil(max(times)/1e6)

    ms_count = 0
    masterlist = []
    for i in range(int(final_ms)):
        ms_bin = times[np.where(np.logical_and(times>i*1e6,times<(i+1)*1e6))]

        if len(ms_bin)<maxcount and len(ms_bin>0):
            masterlist = masterlist + list(ms_bin)
            ms_count+=1

    return np.array(masterlist),ms_count

def find_coincidence(up, down):
    if len(up)>len(down):
        short = down;long = up; sdown = True
    else:
        short = up; long = down; sdown = False
    
    data_list = []
    #print("Analyzing "+str(len(short))+" times...")
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

def write_data(data_arr,ms_vals,filepath):

    filename = filepath.split('/')[-1]


    np.savetxt("/data/p-one/mens/Coincidence_analysis/parsing/coincidences/newJulyData/msCounts/"+filename+"_ms.txt",ms_vals,delimiter="\t")
    np.savetxt("/data/p-one/mens/Coincidence_analysis/parsing/coincidences/newJulyData/coincidenceFiles/"+filename+"_coincs.txt",data_arr,delimiter="\t")



filepath = sys.argv[1]

data = read_in_data(filepath)

uptimes = data[2][np.where(np.logical_and(data[0]==1,data[1]==0))]
downtimes = data[2][np.where(np.logical_and(data[0]==5,data[1]==0))]


parsed_up, up_ms = parse(uptimes)
parsed_down, down_ms = parse(downtimes)

ms_count = [up_ms,down_ms]

coincs = find_coincidence(parsed_up,parsed_down)

write_data(coincs,np.array([[max(ms_count),min(ms_count)]]),filepath)


