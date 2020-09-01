#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3
'''
Functions used for analyzing STRAW data and STRAW GEANT4 simulations

Functions:
    QE(energy)
    read_in_data(filepath)
    parse(times,count_ms=False)
    find_coincidence(up, down)
    find_sim_coincidence(up, down ,urad, drad, uene, dene)
    write_data(data_arr,filepath,filename)
    get_index(array,index,ud)
'''

import sys
sys.path.append('/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/python3.6/site-packages')
import numpy as np
import pandas as pd 
import os
from scipy.interpolate import CubicSpline
import random
import warnings

#get QE function from HamamatsuQE.txt for QE function further down, should work on init
module_dir = os.path.dirname(__file__)
QE_file = np.loadtxt(module_dir+'/HamamatsuQE.txt',delimiter="\t")
QE_func = CubicSpline(QE_file.T[0],QE_file.T[1],extrapolate=False)
del QE_file

#if numpy version > 1.16 need rng for random number generator but as of Aug 2020 Illume uses numpy1.14
try:
    rng = np.random.default_rng()
except AttributeError:
    rng = np.random

#return the QE weight for a given energy
def QE(energy):
    '''
    Returns the Quantum efficiency (QE) value for the input energy

    Parameters:
        energy (array-like): The energy of the detected photon in eV.

    Returns:
        QE (array-like):  The QE of the detected photons for the given energies. If the input energy is not in 
                            the range of 2-4.5 eV, the value will be returned as NaN.
    '''
    return QE_func(energy)


def read_in_data(filepath):
    '''
    Reads a STRAW (or G4 simulation) data file and returns that file as a numpy array

    Parameters:
        filepath (string): The path to the file to be read

    Returns:
        data_array (numpy ndarray): The function will read the file given by the file path assuming it has
                                    been prepared similarly to a STRAW data file (with the delimiter being a '\t' and
                                    no header) and with comments being given with a '#' (this is for functionality with
                                    the simulation files and will not matter for STRAW data unless things change)
    '''
    try:
        data = pd.read_csv(filepath,delimiter='\t',header=None,comment="#")
    except pd.errors.EmptyDataError: exit(1)
    
    return np.array(data)

def parse(times,count_ms=False): #parse by the ms and if we have more than 20 hits in a ms throw out
    '''
    Parse an array of times (in nanoseconds) by the millisecond and return the times that are in low (less than 20/ms)
    millisecond bins

    Parameters:
        times (numpy-array): The list of times from the data (or simulation) file
        
        count_ms (bool): Flag whether to count how many ms bins you are keeping. 
                         Used for getting a frequency of hits/coincidences instead of just a count Default False

    Returns:
        data_list (numpy-array): The list of times where there are less than 20 hits in each time's respective millisecond
        ms_count (integer): If count_ms=True this will return the number of ms bins kept.
    '''
    
    maxcount = 20
 #ms
    final_ms = np.ceil(max(times)/1e6)

    masterlist = []
    ms_count=0
    for i in range(int(final_ms)):
        ms_bin = times[np.where(np.logical_and(times>i*1e6,times<(i+1)*1e6))]

        if len(ms_bin)<maxcount and len(ms_bin)>0:
            masterlist = masterlist + list(ms_bin)
            ms_count+=1
    if count_ms:
        return np.array(masterlist), ms_count
    else:
        return np.array(masterlist)

def find_coincidence(up, down):
    '''
    Find coincidences between two time arrays. Meant for STRAW data files

    Parameters:
        up (numpy array): First array of times. If this is a STRAW data file (or similar) 
                            then this should correspond to the times for the upward facing PMT.

        down (numpy array): Second array of times. If this is a STRAW data file (or similar) then this
                            should correspond to the times for the downward facing PMT.

    Returns:
        data_list (numpy array): nx2 list of time differences for coincident hits with the first column being the
                                    unsorted differences and the second column always being 'up' minus 'down'. 
                                    All time differences will be less than 25ns since this is the cutoff for
                                    consideration as a coincidence
    '''
    if len(up)>len(down):
        short = down;long = up; sdown = True
    else:
        short = up; long = down; sdown = False
    
    data_list = []
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

def find_sim_coincidence(up, down, urad,drad, uene, dene):
    '''
    Find coincidences between two time arrays along with additional information. Meant for GEANT4 simulation files

    Parameters:
        up (numpy-array): The first array of times. If this is a GEANT4 simulation file (or similar) this should
                            correspond to the times for the upward facing PMT.

        down (numpy-array): The second array of times. If this is a GEANT4 simulation file (or similar) this should
                            correspond to the times for the downward facing PMT.

        urad (numpy-array): The list of radial locations for the photons detected at the times listed in the up array.
                            This array should have the same length as the up array.

        drad (numpy-array): The list of radial locations for the photons detected at the times listed in the down array.
                            This array should have the same length as the down array.

        uene (numpy-array): The energies of the photons detected at the times listed in the up array. This array should 
                            have the same length as the up array.

        dene (numpy-array): The energies of the photons detected at the times listed in the down array. This array should
                            have the same length as the down array.


    Returns:
        data_list (numpy-array): An nx6 array of the information on each coincident hit found. The first and second
                                 columns of the array correspond to the time difference between the coincident hit
                                 with the second column always being 'up' minus 'down'. The third and fourth columns 
                                 are the radial locations of the coincident photons detected by the upward and downward facing PMTs,
                                 respectively with the time difference given in the first (or second) column. 
                                 The fifth and sixth columns are the energies of the coincident photons detected in the
                                 upward and downward facing PMTs, respectively with the time difference given in the
                                 first (or second) column. A coincidence is considered as a hit on both the upward and 
                                 downward facing PMTs within a time difference of 25ns.
    '''
    if len(up)>len(down):
        short = down;long = up; sdown = True
        srad = drad; lrad = urad; sene = dene; lene = uene
    else:
        short = up; long = down; sdown = False
        srad = urad; lrad = drad; sene = uene; lene = dene
    
    data_list = []
    for i in range(len(short)):
        t0 = short[i]
        tmin = t0-25;tmax = t0+25 #give 25ns for a 'coincidence'
        potential_list = long[np.logical_and(long>tmin,long<tmax)]
    
        if len(potential_list)>0:
            sradius = srad[i]
            senergy = sene[i]
            for j in range(len(potential_list)):
                diff = t0-potential_list[j]
                lradius = lrad[long==potential_list[j]][0]
                lenergy = lene[long==potential_list[j]][0]
                if (sdown): data_list.append([diff,-diff,lradius,sradius,lenergy,senergy]) #second element is always up-down
                else: data_list.append([diff,diff,sradius,lradius,senergy,lenergy]) #upradius, downradius
    
    data_list = np.array(data_list)
    print(data_list)
    return data_list


def write_data(data_arr,filepath,filename,append=False):
    '''
    Write an array of data out to a file.

    Parameters:
        data_arr (numpy array): The array of data to be written out. Must be 1-D or 2-D.

        filepath (string): The path to the directory you want this data array to be written to.

        filename (string): The name of the file (including extension).

        append (bool): Default False, whether or not you want a new file at filepath with filename or to append
                        to an existing file at filepath

    Returns:
        None, the function will write a file with the given name to the directory specified by the filepath with
        a delimiter of '\t' for 2D arrays.
    '''
    
    if str(filepath[-1])!='/': filepath = str(filepath)+'/'

    if append:
        with open(str(filepath)+str(filename),'a') as f:
            np.savetxt(f,data_arr,delimiter="\t")
    else:
        np.savetxt(str(filepath)+str(filename),data_arr,delimiter="\t")


def get_index(array,index,ud): #ud is a string with 'up' or 'down'
    '''
    Return the values from a specific index of a 2D array for either the upward or downward PMT channels

    Parameters:
        array (numpy-array, STRAW data-like): The array of data which you would like to get the index from. The array
                                                must have the first two columns representing the channel (1 for up, 
                                                5 for down) and the edge (rising=0, falling=1) 

        index (intger): The index for the column of data you would like

        ud (string): Must be either 'up' or 'down'. Denotes which values to return. 'up' will return the values 
                     corresponding to the upward facing pmt (channel 1) and 'down' will return the values corresponding
                     to the downward facing pmt (channel 5). If anything other than 'up' or 'down' is put as this
                     parameter a ValueError will be raised.

    Returns:
        out_array (numpy-array): The column of data specified by the input parameters 
    '''
    
    if str(ud).lower() == "up": 
       ch=1
    elif str(ud).lower() == "down": 
        ch =5
    else:
        raise ValueError("Invalid up/down string entered")
        sys.exit(1)

    return array.T[index][np.logical_and(array.T[0]==ch,array.T[1]==0)]



