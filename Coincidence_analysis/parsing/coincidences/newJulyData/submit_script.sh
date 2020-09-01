universe = vanilla
executable = /cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3
arguments = /data/p-one/mens/Coincidence_analysis/parsing/coincidences/newJulyData/analysis+msCount.py /data/p-one/mens/Coincidence_analysis/parsing/coincidences/newJulyData/newData/file.$(Process) 

output = /data/p-one/mens/condor_logfiles/parsing/newJuly/out/$(Process).out
log = /data/p-one/mens/condor_logfiles/parsing/newJuly/log/$(Process).log
error = /data/p-one/mens/condor_logfiles/parsing/newJuly/err/$(Process).err

request_memory = 2048

periodic_remove = (CommittedTime - CommittedSuspensionTime) > 43200

queue 1350
