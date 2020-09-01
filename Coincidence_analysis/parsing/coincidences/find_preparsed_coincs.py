#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import sys
sys.path.append('/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/python3.6/site-packages')
import numpy as np
import os

def read_in_data(filepath):
    utimes = np.loadtxt(filepath+"_up_parsed.txt",delimiter="\t")
    dtimes = np.loadtxt(filepath+"_down_parsed.txt",delimiter="\t")

    return utimes, dtimes

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
    
    np.savetxt("/data/p-one/mens/Coincidence_analysis/parsing/coincidences/processed/"+filename+"_coincs.txt",data_arr,delimiter="\t")


filepath = sys.argv[1]

up,down = read_in_data(filepath)

coincs = find_coincidence(up,down)
write_data(coincs,filepath)


