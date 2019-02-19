#!/bin/sh

set -x

if [[ ! -d ../exec ]];then
  mkdir ../exec
fi

mac=$(hostname -f)

case $mac in

#---------------------------------------------------------------------------------
# BUILD PROGRAM ON WCOSS Phase 1/2.
#---------------------------------------------------------------------------------

g????.ncep.noaa.gov | t????.ncep.noaa.gov)

  echo NOT ADDED YET
  exit ;;

#---------------------------------------------------------------------------------
# BUILD PROGRAM ON WCOSS CRAY.
#---------------------------------------------------------------------------------

llogin? | slogin?)

  module purge
  module load modules/3.2.6.7
  module load PrgEnv-intel/5.2.56
  module rm intel
  module load intel/16.3.210
  module load cray-mpich/7.2.0
  module load craype-haswell
  module load cray-netcdf
  module load w3nco-intel/2.0.6

  module use /gpfs/hps3/emc/nems/noscrub/emc.nemspara/soft/modulefiles
  module load esmf/7.1.0r

  export FCOMP=ftn
  export FFLAGS="-O0 -g -r8 -i4 -convert big_endian -check bounds -assume byterecl -warn unused"

  make clean
  make
  rc=$?  ;;

#---------------------------------------------------------------------------------
# BUILD PROGRAM ON THEIA.
#---------------------------------------------------------------------------------

tfe??)

  source /apps/lmod/lmod/init/sh
  module purge

  module load intel/15.1.133
  module load impi/5.1.1.109
  module load netcdf/4.3.0

  module use /scratch4/NCEPDEV/nems/noscrub/emc.nemspara/soft/modulefiles
  module load esmf/7.1.0r

  module use -a /scratch3/NCEPDEV/nwprod/lib/modulefiles
  module load w3nco

  export FCOMP=mpiifort
  export FFLAGS="-O0 -g -traceback -r8 -i4 -convert big_endian -check bounds -warn unused -assume byterecl"

  make clean
  make
  rc=$?  ;;

#---------------------------------------------------------------------------------
# BUILD PROGRAM ON DELL.
#---------------------------------------------------------------------------------

m????.ncep.noaa.gov | v????.ncep.noaa.gov )

  module purge
  module use /usrx/local/dev/modulefiles

  module load ips/18.0.1.163
  module load impi/18.0.1
  module load NetCDF/4.5.0
  module load ESMF/7_1_0r
  module load w3nco/2.0.6

  export FCOMP=mpif90
  export FFLAGS="-O0 -g -traceback -r8 -i4 -convert big_endian -check bounds -warn unused -assume byterecl"

  make clean
  make
  rc=$? ;;

esac

exit
