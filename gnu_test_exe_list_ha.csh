#PBS -lselect=1:ncpus=24:model=has      ## Use only one node but access all CPUs available per node.
                                        ## In this case, Haswell was chosen which offers 24 CPUs per node.
#PBS -lwalltime=1:00:00                 ## Keep Wall time under two hours to cut down on queue times.
#PBS -J 1-2:24                          ## Job numbers that will be run in parallel, counting by 24
##PBS -J 1-24:24
set i = $PBS_ARRAY_INDEX
set j = $PBS_ARRAY_INDEX
@ j = $j + 23
echo $i $j
seq $i $j | parallel -j 24 -u --sshloginfile "$PBS_NODEFILE" \
"cd $PWD;./myscript_list_testexe.csh {}"
