#!/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3

import sys
sys.path.append('/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/python3.6/site-packages')
import os
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 

filepath = sys.argv[1]
filename=filepath.split('/')[-1]

d = pd.read_csv(filepath,header=None,delimiter="\t")
channel = np.array(d[0])
edge = np.array(d[1])
time = np.array(d[2])

chlist = [1,5];colorlist = ["C0","C1"]
for i in range(2):
	ch= chlist[i]; col = colorlist[i]
	rise = time[np.where(np.logical_and(channel==ch,edge==0))]
	hist = np.histogram(rise,bins=int(np.ceil(max(time[channel==0])/1e6)))
	plt.close()
	plt.plot(hist[1][:-1]/1e6,hist[0],'x--',color=col,label="channel "+str(ch))
	plt.title(filename)
	plt.xlabel('time [ms]')
	plt.ylabel('Counts')
	plt.legend(loc=0)
	plt.savefig('pdfs/'+filename[:-4]+'_ch'+str(ch)+'.png', format='png', dpi=100)
	plt.ylim(-5,100)
	plt.savefig('ylim_pdfs/'+filename[:-4]+'_ch'+str(ch)+'.png', format='png', dpi=100)


