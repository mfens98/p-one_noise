#!/Users/matthewens/opt/anaconda3/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

readout = pd.read_csv('analysis_QE.txt',delimiter='\t',header=None)
angle = np.array(readout[0])
hits = np.array(readout[1])

plt.figure()
plt.errorbar(angle,hits/hits[0],yerr=np.ones(len(hits))*np.sqrt(hits[0])/hits[0],color="C0",fmt='o--',label="Simulated Acceptance")
plt.plot(angle,np.cos(0.5*angle*np.pi/180)**4,label=r'$\cos\left(0.5\cdot\theta\right)^4$',color="C1")
plt.legend(loc=0)
plt.title("Angular acceptance of simulated sDOM")
plt.xlabel("Angle [deg]")
plt.ylabel("Normalized Intensity")
plt.savefig("acceptanceplot.pdf")
plt.show()
plt.close()
