universe = vanilla
executable = /data/p-one/mens/Coincidence_analysis/parsing/pulse_check/run.sh 

output = /data/p-one/mens/condor_logfiles/pulseCheck1/out/$(Process).out
log = /data/p-one/mens/condor_logfiles/pulseCheck1/log/$(Process).log
error = /data/p-one/mens/condor_logfiles/pulseCheck1/err/$(Process).err

request_memory = 2048

periodic_remove = (CommittedTime - CommittedSuspensionTime) > 43200

queue 
