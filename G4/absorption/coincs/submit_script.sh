universe = vanilla
executable = /data/p-one/mens/G4/absorption/coincs/absCoi.sh
arguments = 8*$(Process) 8*$(Process)+7

output = /data/p-one/mens/condor_logfiles/G4_absorption/out/$(Process).out
log = /data/p-one/mens/condor_logfiles/G4_absorption/log/$(Process).log
error = /data/p-one/mens/condor_logfiles/G4_absorption/err/$(Process).err

request_memory = 1024

periodic_remove = (CommittedTime - CommittedSuspensionTime) > 43200


queue 10000
