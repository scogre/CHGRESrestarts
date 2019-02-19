#!/bin/ksh

#-----------------------------------------------------------
# Run test case on Theia.  MUST BE RUN WITH A 
# MULTIPLE OF SIX MPI TASKS.  Could not get it to
# work otherwise.
#-----------------------------------------------------------

#PBS -l nodes=1:ppn=6
#PBS -l walltime=0:05:00
#PBS -A fv3-cpu
#PBS -q debug
#PBS -N fv3
#PBS -o ./log
#PBS -e ./log

set -x

np=$PBS_NP

source /apps/lmod/lmod/init/ksh
module purge
module load intel/15.1.133
module load impi/5.1.1.109 
module load netcdf/4.3.0

WORKDIR=/scratch3/NCEPDEV/stmp1/$LOGNAME/chgres_fv3
rm -fr $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR

#ln -fs ${PBS_O_WORKDIR}/config.C48.theia.nml ./fort.41
#ln -fs ${PBS_O_WORKDIR}/config.C384.theia.nml ./fort.41
#ln -fs ${PBS_O_WORKDIR}/config.C768.nest.theia.nml ./fort.41
ln -fs ${PBS_O_WORKDIR}/config.C96.nest.theia.nml ./fort.41
#ln -fs ${PBS_O_WORKDIR}/config.C768.stretch.theia.nml ./fort.41
#ln -fs ${PBS_O_WORKDIR}/config.C1152.theia.nml ./fort.41

mpirun -np $np ${PBS_O_WORKDIR}/../exec/global_chgres.exe

exit 0
