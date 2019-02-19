 module program_setup

 implicit none

 private

 character(len=500), public      :: data_dir_input_grid = "NULL"
 character(len=500), public      :: fix_dir_target_grid = "NULL"
 character(len=500), public      :: mosaic_file_input_grid = "NULL"
 character(len=500), public      :: mosaic_file_target_grid = "NULL"
 character(len=500), public      :: orog_dir_input_grid = "NULL"
 character(len=500), public      :: orog_files_input_grid(6) = "NULL"
 character(len=500), public      :: orog_dir_target_grid = "NULL"
 character(len=500), public      :: orog_files_target_grid(6) = "NULL"
 character(len=500), public      :: vcoord_file_target_grid = "NULL"
 character(len=6),   public      :: cres_target_grid = "      "

 integer, public                 :: cycle_mon = -999
 integer, public                 :: cycle_day = -999
 integer, public                 :: cycle_hour = -999
 integer, public                 :: regional = 0
 integer, public                 :: halo = 0
 integer, public                 :: halo_p1 = 0

 logical, public                 :: convert_atm = .false.
 logical, public                 :: convert_nst = .false.
 logical, public                 :: convert_sfc = .false.
 logical, public                 :: gfdl_mp = .false.
 logical, public                 :: restart_file = .false.

 real, allocatable, public       :: drysmc_input(:), drysmc_target(:)
 real, allocatable, public       :: maxsmc_input(:), maxsmc_target(:)
 real, allocatable, public       :: refsmc_input(:), refsmc_target(:)
 real, allocatable, public       :: wltsmc_input(:), wltsmc_target(:)
 real, allocatable, public       :: bb_target(:),    satpsi_target(:)

 public :: read_setup_namelist
 public :: calc_soil_params_driver

 contains

 subroutine read_setup_namelist

 implicit none

 integer                     :: is, ie, ierr

 include 'mpif.h'

 namelist /config/ mosaic_file_target_grid, &
                   fix_dir_target_grid,     &
                   orog_dir_target_grid,    &
                   orog_files_target_grid,  &
                   mosaic_file_input_grid,  &
                   orog_dir_input_grid,     &
                   orog_files_input_grid,   &
                   data_dir_input_grid,     &
                   vcoord_file_target_grid, &
                   cycle_mon, cycle_day,    &
                   cycle_hour, convert_atm, &
                   convert_nst, convert_sfc, &
                   regional, restart_file

 print*,"- READ SETUP NAMELIST"

 open(41, file="./fort.41", iostat=ierr, err=900)
 read(41, nml=config, iostat=ierr, err=901)
 close (41)

 orog_dir_target_grid = trim(orog_dir_target_grid) // '/'
 orog_dir_input_grid = trim(orog_dir_input_grid) // '/'

 is = index(mosaic_file_target_grid, "/", .true.)
 ie = index(mosaic_file_target_grid, "_mosaic")

 if (is == 0 .or. ie == 0) then
   print*,'bad cres '
   call mpi_abort
 endif
   
 cres_target_grid = mosaic_file_target_grid(is+1:ie-1)

!-------------------------------------------------------------------------
! Flag for processing stand-alone regional grid.  When '1', 
! remove halo from atmospheric and surface data and output
! atmospheric lateral boundary condition file. When '2',
! create lateral boundary file only.  When '0' (the default),
! process normally as a global grid.
!-------------------------------------------------------------------------

 if (regional > 0) then
   halo = 3
   halo_p1 = halo + 1
   print*,"- PROCESSING A REGIONAL NEST WITH A HALO OF ",halo
 endif

 return

 900 print*,'- FATAL ERROR OPENING CONFIG NAMELIST'
 print*,'- IOSTAT IS: ', ierr
 call mpi_abort(mpi_comm_world, 11, ierr)

 901 print*,'- FATAL ERROR READING CONFIG NAMELIST'
 print*,'- IOSTAT IS: ', ierr
 call mpi_abort(mpi_comm_world, 11, ierr)

 end subroutine read_setup_namelist

 subroutine calc_soil_params_driver(localpet)

 implicit none

 integer, intent(in)       :: localpet

 integer, parameter        :: num_statsgo = 16
 real, parameter           :: smlow_statsgo = 0.5
 real, parameter           :: smhigh_statsgo = 6.0

 integer                   :: num_soil_cats

 real                      :: bb_statsgo(num_statsgo)
 real                      :: maxsmc_statsgo(num_statsgo)
 real                      :: satdk_statsgo(num_statsgo)
 real                      :: satpsi_statsgo(num_statsgo)

 real, allocatable         :: bb(:)
 real                      :: smlow, smhigh
 real, allocatable         :: f11(:)
 real, allocatable         :: satdk(:)
 real, allocatable         :: satpsi(:)
 real, allocatable         :: satdw(:)

 data bb_statsgo /4.05, 4.26, 4.74, 5.33, 5.33, 5.25, &
            6.77, 8.72, 8.17, 10.73, 10.39, 11.55, &
            5.25, -9.99, 4.05, 4.26/

 data maxsmc_statsgo /0.395, 0.421, 0.434, 0.476, 0.476, 0.439, &
              0.404, 0.464, 0.465, 0.406, 0.468, 0.457, &
              0.464, -9.99, 0.200, 0.421/

 data satdk_statsgo /1.7600e-4, 1.4078e-5, 5.2304e-6, 2.8089e-6, 2.8089e-6, &
             3.3770e-6, 4.4518e-6, 2.0348e-6, 2.4464e-6, 7.2199e-6, &
             1.3444e-6, 9.7384e-7, 3.3770e-6,     -9.99, 1.4078e-5, &
             1.4078e-5/

 data satpsi_statsgo /0.0350, 0.0363, 0.1413, 0.7586, 0.7586, 0.3548, &
              0.1349, 0.6166, 0.2630, 0.0977, 0.3236, 0.4677, &
              0.3548, -9.99,  0.0350, 0.0363/

! input grid

 num_soil_cats = num_statsgo

 allocate(maxsmc_input(num_soil_cats))
 allocate(wltsmc_input(num_soil_cats))
 allocate(drysmc_input(num_soil_cats))
 allocate(refsmc_input(num_soil_cats))
 allocate(bb(num_soil_cats))
 allocate(satdk(num_soil_cats))
 allocate(satpsi(num_soil_cats))
 allocate(satdw(num_soil_cats))
 allocate(f11(num_soil_cats))

 smlow  = smlow_statsgo
 smhigh = smhigh_statsgo
 maxsmc_input = maxsmc_statsgo
 bb     = bb_statsgo
 satdk  = satdk_statsgo
 satpsi = satpsi_statsgo

 call calc_soil_params(num_soil_cats, smlow, smhigh, satdk, maxsmc_input, &
                       bb, satpsi, satdw, f11, refsmc_input, drysmc_input, wltsmc_input)

 deallocate(bb, satdk, satpsi, satdw, f11)

 if (localpet == 0) print*,'maxsmc input grid ',maxsmc_input

! target grid

 num_soil_cats = num_statsgo

 allocate(maxsmc_target(num_soil_cats))
 allocate(wltsmc_target(num_soil_cats))
 allocate(drysmc_target(num_soil_cats))
 allocate(refsmc_target(num_soil_cats))
 allocate(bb_target(num_soil_cats))
 allocate(satpsi_target(num_soil_cats))
 allocate(satdk(num_soil_cats))
 allocate(satdw(num_soil_cats))
 allocate(f11(num_soil_cats))

 smlow  = smlow_statsgo
 smhigh = smhigh_statsgo
 maxsmc_target = maxsmc_statsgo
 bb_target     = bb_statsgo
 satdk  = satdk_statsgo
 satpsi_target = satpsi_statsgo

 call calc_soil_params(num_soil_cats, smlow, smhigh, satdk, maxsmc_target, &
                       bb_target, satpsi_target, satdw, f11, refsmc_target, drysmc_target, wltsmc_target)

 deallocate(satdk, satdw, f11)

 if (localpet == 0) print*,'maxsmc target grid ',maxsmc_target

 end subroutine calc_soil_params_driver

 subroutine calc_soil_params(num_soil_cats, smlow, smhigh, satdk,  &
            maxsmc, bb, satpsi, satdw, f11, refsmc, drysmc, wltsmc)

 implicit none

 integer, intent(in)            :: num_soil_cats

 real, intent(in)               :: smlow, smhigh
 real, intent(in)               :: bb(num_soil_cats)
 real, intent(in)               :: maxsmc(num_soil_cats)
 real, intent(in)               :: satdk(num_soil_cats)
 real, intent(in)               :: satpsi(num_soil_cats)

 real, intent(out)              :: f11(num_soil_cats)
 real, intent(out)              :: satdw(num_soil_cats)
 real, intent(out)              :: refsmc(num_soil_cats)
 real, intent(out)              :: drysmc(num_soil_cats)
 real, intent(out)              :: wltsmc(num_soil_cats)

 integer                        :: i

 real                      :: refsmc1
 real                      :: wltsmc1

 satdw = 0.0
 f11   = 0.0
 refsmc = 0.0
 wltsmc = 0.0
 drysmc = 0.0

 do i = 1, num_soil_cats

   if (maxsmc(i) > 0.0) then

   SATDW(I)  = BB(I)*SATDK(I)*(SATPSI(I)/MAXSMC(I))
   F11(I) = ALOG10(SATPSI(I)) + BB(I)*ALOG10(MAXSMC(I)) + 2.0
   REFSMC1 = MAXSMC(I)*(5.79E-9/SATDK(I)) **(1.0/(2.0*BB(I)+3.0))
   REFSMC(I) = REFSMC1 + (MAXSMC(I)-REFSMC1) / SMHIGH
   WLTSMC1 = MAXSMC(I) * (200.0/SATPSI(I))**(-1.0/BB(I))
   WLTSMC(I) = WLTSMC1 - SMLOW * WLTSMC1

!----------------------------------------------------------------------
!  CURRENT VERSION DRYSMC VALUES THAT EQUATE TO WLTSMC.
!  FUTURE VERSION COULD LET DRYSMC BE INDEPENDENTLY SET VIA NAMELIST.
!----------------------------------------------------------------------

   DRYSMC(I) = WLTSMC(I)

   end if

 END DO

 end subroutine calc_soil_params

 end module program_setup
