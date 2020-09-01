#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import os
import re

codir='/data/p-one/mens/Coincidence_analysis/parsing/pulse_check/coincidenceFiles/'

files=os.listdir(codir)

with open(codir+'../checkFiles.txt') as f:
	datafiles=f.readlines()
datafiles=[x.strip() for x in datafiles]
mastercount=20
for f in files:
	mastercount+=1
	print("\n\nStarting file %s, number %d" % (f,mastercount))
	name=f.split('.')[0]+'.txt'



	coincs = np.loadtxt(codir+f,delimiter="\t")

	dfile = [i for i in datafiles if name in i] #get datafile path
	try:
		d=pd.read_csv(dfile[0],delimiter="\t",header=None)
	except IndexError as exp:
		print(exp)
		print(dfile)
		continue
	ad = np.array(d)
	data=ad[ad[:,2].argsort()] #get array sorted by times

	count=0
	for c in coincs:
		print("Looking at coinc number %d " % count)
		upulse = []
		dpulse = []
		for pulse in [c[0],c[1]]:
			start = np.where(data.T[2]==pulse)[0][0]
			i = start; flag=True; start_point = data[start] 
			while flag: 
				if (data[i][0]==start_point[0] and data[i][1]==1): flag=False

				if (data[i][0]-start_point[0])<=3 and (data[i][0]-start_point[0])>=0:
					if start_point[0] == 1: upulse.append(data[i])
					elif start_point[0]==5: dpulse.append(data[i])
				i+=1
		upulse=np.array(upulse)
		dpulse=np.array(dpulse)
		#check ordering and just save times
		try:
			if np.array_equal(np.sort(upulse.T[0][upulse.T[1]==0]),upulse.T[0][upulse.T[1]==0]) and \
			   np.array_equal(np.sort(upulse.T[0][upulse.T[1]==1])[::-1],upulse.T[0][upulse.T[1]==1]):

			   utimes=upulse.T[2]
			else:
				print("Uptimes failed")
				print(upulse)
		except IndexError:
			print("Uptimes Index Error")
			print(upulse)
			count+=1
			continue

		try:
			if np.array_equal(np.sort(dpulse.T[0][dpulse.T[1]==0]),dpulse.T[0][dpulse.T[1]==0]) and \
		       np.array_equal(np.sort(dpulse.T[0][dpulse.T[1]==1])[::-1],dpulse.T[0][dpulse.T[1]==1]):

		   		dtimes=dpulse.T[2]
			else:
				print("Downtimes failed")
				print(dpulse)
		except IndexError:
			print("Downtimes Index Error")
			print(dpulse)
			count+=1
			continue

		plt.close()
		strike=0
		try:
			if len(utimes)/2 % 1 !=0 : raise TypeError
			plt.plot(utimes,np.array(list(range(int(len(utimes)/2)))+list(range(int(len(utimes)/2)))[::-1])+1,'o--',label="Upward PMT",color="C0")
		except NameError: #utimes was not defined since test failed
			strike+=1; pass
		except TypeError: #length of utimes was not even so there was an unread value
			strike+=1; print("Uptime length %d\n" % len(utimes),upulse); pass
		try:
			if len(dtimes)/2 % 1 != 0 : raise TypeError
			plt.plot(dtimes,np.array(list(range(int(len(dtimes)/2)))+list(range(int(len(dtimes)/2)))[::-1])+1,'o--',label="Downward PMT",color="C1")
		except NameError: #dtimes was not defined since test failed
			strike+=1; pass
		except TypeError: #length of dtimes was not even so there was an unread value
			strike+=1; print("Downtime length %d\n" % len(dtimes),dpulse); pass
		if strike==2: print("Too many strikes"); count+=1; continue
		plt.xlabel("Time [ns]")
		plt.ylabel("Threshold")
		plt.title(name+'_'+str(count))
		plt.legend(loc=0)
		plt.savefig(codir+'../pngs/'+name+'_'+str(count)+'.png',format='png',dpi=100)

		count+=1
