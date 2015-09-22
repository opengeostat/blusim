#!/bin/bash 
# The total number of processes:
# 64 MPI processes * 2 OpenMP threads = 128 cores
#BSUB -n 16
#BSUB -oo job_%J.out
#BSUB -eo job_%J.err
# It will allocate 1 mpi proc per node:
#BSUB -R"span[ptile=16]"
#BSUB -x
#BSUB -J parblusim
#BSUB -W 00:30
#BSUB -q bsc_debug

# You can choose the parallel environment through
# modules
#module load intel
#module load openmpi


export OMP_NUM_THREADS=1

SIZE=small
ID=`date +%s`

echo '--| LSF: JOB SUBMITTED WITH SRUN AT: ' `date`
#/usr/bin/time mpirun /home/bsc21/bsc21021/Alya_repsol/Alya_MN/Executables/unix/Alya.x Invo4 > /gpfs/scratch/bsc21/bsc21021/outputs/output-MN-${ID}.txt
#/usr/bin/time /home/bsc21/bsc21021/parblusim/blusim.exe blusim-large.par > /home/bsc21/bsc21021/parblusim/output-MN-${ID}-${OMP_NUM_THREADS}.txt
for i in 1 2 3 4 5
do
/usr/bin/time /home/bsc21/bsc21021/parblusim/blusim.exe blusim-${SIZE}.par > /home/bsc21/bsc21021/parblusim/output-MN-${ID}-${OMP_NUM_THREADS}.txt
done
#mpirun -v /gpfs/apps/MN3/VALGRIND/3.8.1/bin/valgrind --error-limit=no --leak-check=full /home/bsc21/bsc21021/Alya_merge/Executables/unix/Alya.g Invo4
echo '--| LSF: JOB FINISHED AT: ' `date`

#./moveSolution.sh ${OMP_NUM_THREADS} ${SIZE}
