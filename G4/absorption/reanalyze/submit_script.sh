universe = vanilla
executable = /data/p-one/mens/G4/absorption/reanalyze/run.sh 
arguments =  7767*$(Process) 7767*$(Process)+7766

output = /data/p-one/mens/condor_logfiles/reanalysis2/out/$(Process).out
log = /data/p-one/mens/condor_logfiles/reanalysis2/log/$(Process).log
error = /data/p-one/mens/condor_logfiles/reanalysis2/err/$(Process).err

request_memory = 1024

periodic_remove = (CommittedTime - CommittedSuspensionTime) > 43200


queue 103
