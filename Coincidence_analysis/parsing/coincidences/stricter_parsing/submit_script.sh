universe = vanilla
executable = /cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/python3
arguments = /data/p-one/mens/Coincidence_analysis/parsing/coincidences/stricter/parsed_analysis.py /data/p-one/mens/Coincidence_analysis/parsing/datafiles/file.$(Process) 

output = /data/p-one/mens/condor_logfiles/parsing/stricter/out/$(Process).out
log = /data/p-one/mens/condor_logfiles/parsing/stricter/log/$(Process).log
error = /data/p-one/mens/condor_logfiles/parsing/stricter/err/$(Process).err

request_memory = 2048

periodic_remove = (CommittedTime - CommittedSuspensionTime) > 43200

queue 2395
