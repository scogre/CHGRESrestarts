 module atmosphere

 use esmf

 use input_data, only                : lev_input, &
                                       levp_input, &
                                       ntracer_input, &
                                       liq_wat_input_grid, &
                                       o3mr_input_grid, &
                                       rwmr_input_grid, &
                                       snmr_input_grid, &
                                       icmr_input_grid, &
                                       grle_input_grid, &
                                       clda_input_grid, &
                                       dzdt_input_grid, &
                                       ps_input_grid, &
                                       spec_hum_input_grid, &
                                       wind_input_grid,   &
                                       temp_input_grid,   &
                                       pres_input_grid,   &
                                       read_input_atm_data, &
                                       cleanup_input_atm_data

 use model_grid, only                : target_grid,  &
                                       latitude_s_target_grid,  &
                                       longitude_s_target_grid, &
                                       latitude_w_target_grid,  &
                                       longitude_w_target_grid, &
                                       terrain_input_grid,      &
                                       terrain_target_grid

 use program_setup, only             : vcoord_file_target_grid, &
                                       regional, gfdl_mp

 implicit none

 private

 integer, public                    :: lev_target
 integer, public                    :: levp_target
 integer, public                    :: ntracer_target
 integer, public                    :: nvcoord_target

 real, allocatable, public          :: vcoord_target(:,:)

 type(esmf_field), public           :: delp_target_grid
 type(esmf_field), public           :: dzdt_target_grid
 type(esmf_field)                   :: dzdt_b4adj_target_grid
 type(esmf_field), public           :: liq_wat_target_grid
 type(esmf_field)                   :: liq_wat_b4adj_target_grid
 type(esmf_field), public           :: o3mr_target_grid
 type(esmf_field)                   :: o3mr_b4adj_target_grid
 type(esmf_field)                   :: ps_b4adj_target_grid
 type(esmf_field), public           :: ps_target_grid
 type(esmf_field)                   :: pres_target_grid
 type(esmf_field)                   :: pres_b4adj_target_grid
 type(esmf_field), public           :: spec_hum_target_grid
 type(esmf_field)                   :: spec_hum_b4adj_target_grid
 type(esmf_field), public           :: grle_target_grid
 type(esmf_field)                   :: grle_b4adj_target_grid
 type(esmf_field), public           :: clda_target_grid
 type(esmf_field)                   :: clda_b4adj_target_grid
 type(esmf_field), public           :: icmr_target_grid
 type(esmf_field)                   :: icmr_b4adj_target_grid
 type(esmf_field), public           :: rwmr_target_grid
 type(esmf_field)                   :: rwmr_b4adj_target_grid
 type(esmf_field), public           :: snmr_target_grid
 type(esmf_field)                   :: snmr_b4adj_target_grid
 type(esmf_field), public           :: temp_target_grid
 type(esmf_field)                   :: temp_b4adj_target_grid
 type(esmf_field)                   :: terrain_interp_to_target_grid
 type(esmf_field), public           :: u_s_target_grid
 type(esmf_field), public           :: v_s_target_grid
 type(esmf_field)                   :: wind_b4adj_target_grid
 type(esmf_field)                   :: wind_target_grid
 type(esmf_field)                   :: wind_s_target_grid
 type(esmf_field), public           :: u_w_target_grid
 type(esmf_field), public           :: v_w_target_grid
 type(esmf_field)                   :: wind_w_target_grid
 type(esmf_field), public           :: zh_target_grid

 public :: atmosphere_driver

 contains

!-----------------------------------------------------------------------------------
! Driver routine for atmospheric fields.
!-----------------------------------------------------------------------------------

 subroutine atmosphere_driver(localpet)

 implicit none

 include 'mpif.h'

 integer, intent(in)                :: localpet

 integer                            :: isrctermprocessing
 integer                            :: rc

 type(esmf_regridmethod_flag)       :: method
 type(esmf_routehandle)             :: regrid_bl

!-----------------------------------------------------------------------------------
! Read atmospheric fields on the input grid.
!-----------------------------------------------------------------------------------

 call read_input_atm_data(localpet)

!-----------------------------------------------------------------------------------
! Initially, don't change number of tracers.
!-----------------------------------------------------------------------------------
 
 ntracer_target = ntracer_input

!-----------------------------------------------------------------------------------
! Read vertical coordinate info for target grid.
!-----------------------------------------------------------------------------------

 call read_vcoord_info

!-----------------------------------------------------------------------------------
! Create target grid field objects to hold data before vertical adjustment.
!-----------------------------------------------------------------------------------

 call create_atm_b4adj_esmf_fields

!-----------------------------------------------------------------------------------
! Horizontally interpolate.
!-----------------------------------------------------------------------------------

 method=ESMF_REGRIDMETHOD_BILINEAR

 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore FOR ATMOSPHERIC FIELDS."
 call ESMF_FieldRegridStore(spec_hum_input_grid, &
                            spec_hum_b4adj_target_grid, &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_bl, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR TEMPERATURE."
 call ESMF_FieldRegrid(temp_input_grid, &
                       temp_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR PRESSURE."
 call ESMF_FieldRegrid(pres_input_grid, &
                       pres_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR SPECIFIC HUMIDITY."
 call ESMF_FieldRegrid(spec_hum_input_grid, &
                       spec_hum_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR VERTICAL VELOCITY."
 call ESMF_FieldRegrid(dzdt_input_grid, &
                       dzdt_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR OZONE MIXING RATIO."
 call ESMF_FieldRegrid(o3mr_input_grid, &
                       o3mr_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR LIQUID WATER MIXING RATIO."
 call ESMF_FieldRegrid(liq_wat_input_grid, &
                       liq_wat_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort




 if (gfdl_mp) then

   print*,"- CALL Field_Regrid FOR RAIN WATER MIXING RATIO."
   call ESMF_FieldRegrid(rwmr_input_grid, &
                         rwmr_b4adj_target_grid, &
                         routehandle=regrid_bl, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL Field_Regrid FOR SNOW WATER MIXING RATIO."
   call ESMF_FieldRegrid(snmr_input_grid, &
                         snmr_b4adj_target_grid, &
                         routehandle=regrid_bl, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL Field_Regrid FOR ICE WATER MIXING RATIO."
   call ESMF_FieldRegrid(icmr_input_grid, &
                         icmr_b4adj_target_grid, &
                         routehandle=regrid_bl, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL Field_Regrid FOR GRAUPEL."
   call ESMF_FieldRegrid(grle_input_grid, &
                         grle_b4adj_target_grid, &
                         routehandle=regrid_bl, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL Field_Regrid FOR CLOUD AMT"
   call ESMF_FieldRegrid(clda_input_grid, &
                         clda_b4adj_target_grid, &
                         routehandle=regrid_bl, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort


 endif

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore FOR SURFACE PRESSURE."
 call ESMF_FieldRegridStore(ps_input_grid, &
                            ps_b4adj_target_grid, &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_bl, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR SURFACE PRESSURE."
 call ESMF_FieldRegrid(ps_input_grid, &
                       ps_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR TERRAIN."
 call ESMF_FieldRegrid(terrain_input_grid, &
                       terrain_interp_to_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore FOR 3D-WIND."
 call ESMF_FieldRegridStore(wind_input_grid, &
                            wind_b4adj_target_grid, &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_bl, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR 3-D WIND."
 call ESMF_FieldRegrid(wind_input_grid, &
                       wind_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

!-----------------------------------------------------------------------------------
! Deallocate input fields.
!-----------------------------------------------------------------------------------

 call cleanup_input_atm_data

!-----------------------------------------------------------------------------------
! Create target grid field objects to hold data after vertical interpolation.
!-----------------------------------------------------------------------------------

 call create_atm_esmf_fields

!-----------------------------------------------------------------------------------
! Adjust surface pressure for terrain differences.
!-----------------------------------------------------------------------------------

 call newps(localpet)

!-----------------------------------------------------------------------------------
! Compute 3-d pressure based on adjusted surface pressure.
!-----------------------------------------------------------------------------------

 call newpr1(localpet)

!-----------------------------------------------------------------------------------
! Vertically interpolate.
!-----------------------------------------------------------------------------------

 call vintg(localpet)

!-----------------------------------------------------------------------------------
! Compute height.
!-----------------------------------------------------------------------------------

 call compute_zh

!-----------------------------------------------------------------------------------
! Free up memory.
!-----------------------------------------------------------------------------------

 call cleanup_target_atm_b4adj_data

!-----------------------------------------------------------------------------------
! Interpolate winds to 'd' grid.
!-----------------------------------------------------------------------------------

 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore FOR 3D-WIND WEST EDGE."
 call ESMF_FieldRegridStore(wind_target_grid, &
                            wind_w_target_grid, &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
                            routehandle=regrid_bl, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR 3-D WIND WEST EDGE."
 call ESMF_FieldRegrid(wind_target_grid, &
                       wind_w_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore FOR 3D-WIND SOUTH EDGE."
 call ESMF_FieldRegridStore(wind_target_grid, &
                            wind_s_target_grid, &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
                            routehandle=regrid_bl, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL Field_Regrid FOR 3-D WIND SOUTH EDGE."
 call ESMF_FieldRegrid(wind_target_grid, &
                       wind_s_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

!-----------------------------------------------------------------------------------
! Convert from 3-d to 2-d cartesian winds.
!-----------------------------------------------------------------------------------

 call convert_winds

!-----------------------------------------------------------------------------------
! Write target data to file.
!-----------------------------------------------------------------------------------

 call write_fv3_atm_header_netcdf(localpet)
 if (regional <= 1) call write_fv3_atm_data_netcdf(localpet)
 if (regional >= 1) call write_fv3_atm_bndy_data_netcdf(localpet)

!-----------------------------------------------------------------------------------
! Free up memory.
!-----------------------------------------------------------------------------------

 call cleanup_target_atm_data

 end subroutine atmosphere_driver

!-----------------------------------------------------------------------------------
! Create target grid field objects to hold data before vertical interpolation.
! These will be defined with the same number of vertical levels as the input grid.
!-----------------------------------------------------------------------------------

 subroutine create_atm_b4adj_esmf_fields

 implicit none

 integer                          :: rc

 print*,"- CALL FieldCreate FOR TARGET GRID SPEC HUMIDITY BEFORE ADJUSTMENT."
 spec_hum_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID TEMPERATURE BEFORE ADJUSTMENT."
 temp_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID PRESSURE BEFORE ADJUSTMENT."
 pres_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID OZONE BEFORE ADJUSTMENT."
 o3mr_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID LIQ WATER BEFORE ADJUSTMENT."
 liq_wat_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID VERTICAL VELOCITY BEFORE ADJUSTMENT."
 dzdt_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID UNSTAGGERED WINDS BEFORE ADJUSTMENT."
 wind_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1,1/), &
                                   ungriddedUBound=(/lev_input,3/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET TERRAIN."
 terrain_interp_to_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET SURFACE PRESSURE BEFORE ADJUSTMENT."
 ps_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort






 if (gfdl_mp) then

   print*,"- CALL FieldCreate FOR TARGET GRAUPEL BEFORE ADJUSTMENT."
   grle_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL FieldCreate FOR TARGET CLOUD AMT BEFORE ADJUSTMENT."
   clda_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL FieldCreate FOR TARGET ICE WATER BEFORE ADJUSTMENT."
   icmr_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL FieldCreate FOR TARGET RAIN WATER BEFORE ADJUSTMENT."
   rwmr_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL FieldCreate FOR TARGET SNOW WATER BEFORE ADJUSTMENT."
   snmr_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 end if

 end subroutine create_atm_b4adj_esmf_fields

!-----------------------------------------------------------------------------------
! Create target grid field objects.
!-----------------------------------------------------------------------------------

 subroutine create_atm_esmf_fields

 implicit none

 integer                          :: rc

 print*,"- CALL FieldCreate FOR TARGET GRID TEMPERATURE."
 temp_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID PRESSURE."
 pres_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID VERTICAL VELOCITY."
 dzdt_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID DELP."
 delp_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID SPECIFIC HUMIDITY."
 spec_hum_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET GRID OZONE MIXING RATIO."
 o3mr_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET HEIGHT."
 zh_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/levp_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET LIQUID WATER MIXING RATIO."
 liq_wat_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET UNSTAGGERED 3D-WIND."
 wind_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1,1/), &
                                   ungriddedUBound=(/lev_target,3/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET U_S."
 u_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET V_S."
 v_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET 3D-WIND_S."
 wind_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   ungriddedLBound=(/1,1/), &
                                   ungriddedUBound=(/lev_target,3/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET U_W."
 u_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET V_W."
 v_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET 3D-WIND_W."
 wind_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   ungriddedLBound=(/1,1/), &
                                   ungriddedUBound=(/lev_target,3/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 print*,"- CALL FieldCreate FOR TARGET SURFACE PRESSURE."
 ps_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort





 if (gfdl_mp) then

   print*,"- CALL FieldCreate FOR TARGET GRID RAIN WATER MIXING RATIO."
   rwmr_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL FieldCreate FOR TARGET GRID SNOW WATER MIXING RATIO."
   snmr_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL FieldCreate FOR TARGET GRID ICE WATER MIXING RATIO."
   icmr_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL FieldCreate FOR TARGET GRID GRAUPEL RATIO."
   grle_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

   print*,"- CALL FieldCreate FOR TARGET GRID CLOUD AMT."
   clda_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 endif

 end subroutine create_atm_esmf_fields

 subroutine convert_winds
 
 implicit none

 integer                         :: clb(4), cub(4)
 integer                         :: i, j, k, rc

 real(esmf_kind_r8), pointer     :: latptr(:,:)
 real(esmf_kind_r8), pointer     :: lonptr(:,:)
 real(esmf_kind_r8), pointer     :: uptr(:,:,:)
 real(esmf_kind_r8), pointer     :: vptr(:,:,:)
 real(esmf_kind_r8), pointer     :: windptr(:,:,:,:)
 real(esmf_kind_r8)              :: latrad, lonrad

!-----------------------------------------------------------------------------------
! Convert from 3-d cartesian to 2-cartesian winds
!-----------------------------------------------------------------------------------

 print*,'top of convert winds'

 print*,"- CALL FieldGet FOR 3-D WIND_S."
 call ESMF_FieldGet(wind_s_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=windptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 print*,"- CALL FieldGet FOR U_S."
 call ESMF_FieldGet(u_s_target_grid, &
                    farrayPtr=uptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 print*,"- CALL FieldGet FOR V_S."
 call ESMF_FieldGet(v_s_target_grid, &
                    farrayPtr=vptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 print*,"- CALL FieldGet FOR LATITUDE_S."
 call ESMF_FieldGet(latitude_s_target_grid, &
                    farrayPtr=latptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 print*,"- CALL FieldGet FOR LONGITUDE_S."
 call ESMF_FieldGet(longitude_s_target_grid, &
                    farrayPtr=lonptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     latrad = latptr(i,j) * acos(-1.) / 180.0
     lonrad = lonptr(i,j) * acos(-1.) / 180.0
     do k = clb(3), cub(3)
       uptr(i,j,k) = windptr(i,j,k,1) * cos(lonrad) + windptr(i,j,k,2) * sin(lonrad)
       vptr(i,j,k) = -windptr(i,j,k,1) * sin(latrad) * sin(lonrad) + &
                      windptr(i,j,k,2) * sin(latrad) * cos(lonrad) + &
                      windptr(i,j,k,3) * cos(latrad)
     enddo
   enddo
 enddo

 print*,"- CALL FieldGet FOR 3-D WIND_W."
 call ESMF_FieldGet(wind_w_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=windptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 print*,"- CALL FieldGet FOR U_W."
 call ESMF_FieldGet(u_w_target_grid, &
                    farrayPtr=uptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 print*,"- CALL FieldGet FOR V_W."
 call ESMF_FieldGet(v_w_target_grid, &
                    farrayPtr=vptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 print*,"- CALL FieldGet FOR LATITUDE_W."
 call ESMF_FieldGet(latitude_w_target_grid, &
                    farrayPtr=latptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 print*,"- CALL FieldGet FOR LONGITUDE_W."
 call ESMF_FieldGet(longitude_w_target_grid, &
                    farrayPtr=lonptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     latrad = latptr(i,j) * acos(-1.) / 180.0
     lonrad = lonptr(i,j) * acos(-1.) / 180.0
     do k = clb(3), cub(3)
       uptr(i,j,k) = windptr(i,j,k,1) * cos(lonrad) + windptr(i,j,k,2) * sin(lonrad)
       vptr(i,j,k) = -windptr(i,j,k,1) * sin(latrad) * sin(lonrad) + &
                      windptr(i,j,k,2) * sin(latrad) * cos(lonrad) + &
                      windptr(i,j,k,3) * cos(latrad)
     enddo
   enddo
 enddo

 end subroutine convert_winds

 subroutine newpr1(localpet)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK                                    
!                                                                       
! SUBPROGRAM:    NEWPR1      COMPUTE MODEL PRESSURES                    
!   PRGMMR: JUANG          ORG: W/NMC23     DATE: 2005-04-11            
!   PRGMMR: Fanglin Yang   ORG: W/NMC23     DATE: 2006-11-28            
!   PRGMMR: S. Moorthi     ORG: NCEP/EMC    DATE: 2006-12-12            
!   PRGMMR: S. Moorthi     ORG: NCEP/EMC    DATE: 2007-01-02            
!                                                                       
! ABSTRACT: COMPUTE MODEL PRESSURES.                                    
!                                                                       
! PROGRAM HISTORY LOG:                                                  
! 2005-04-11  HANN_MING HENRY JUANG    hybrid sigma, sigma-p, and sigma-
!                                                                       
! USAGE:    CALL NEWPR1(IM,IX,KM,KMP,IDVC,IDSL,NVCOORD,VCOORD,PP,TP,QP,P
!   INPUT ARGUMENT LIST:                                                
!     IM           INTEGER NUMBER OF POINTS TO COMPUTE                  
!     KM           INTEGER NUMBER OF LEVELS                             
!     IDVC         INTEGER VERTICAL COORDINATE ID                       
!                  (1 FOR SIGMA AND 2 FOR HYBRID)                       
!     IDSL         INTEGER TYPE OF SIGMA STRUCTURE                      
!                  (1 FOR PHILLIPS OR 2 FOR MEAN)                       
!     NVCOORD      INTEGER NUMBER OF VERTICAL COORDINATES               
!     VCOORD       REAL (KM+1,NVCOORD) VERTICAL COORDINATE VALUES       
!                  FOR IDVC=1, NVCOORD=1: SIGMA INTERFACE               
!                  FOR IDVC=2, NVCOORD=2: HYBRID INTERFACE A AND B      
!                  FOR IDVC=3, NVCOORD=3: JUANG GENERAL HYBRID INTERFACE
!                     AK  REAL (KM+1) HYBRID INTERFACE A                
!                     BK  REAL (KM+1) HYBRID INTERFACE B                
!     PS           REAL (IX) SURFACE PRESSURE (PA)                      
!   OUTPUT ARGUMENT LIST:                                               
!     PM           REAL (IX,KM) MID-LAYER PRESSURE (PA)                 
!     DP           REAL (IX,KM) LAYER DELTA PRESSURE (PA)               
!                                                                       
! ATTRIBUTES:                                                           
!   LANGUAGE: FORTRAN                                                   
!                                                                       
!C$$$                                                                   
 implicit none 

 integer, intent(in) :: localpet

 integer                         :: idsl, idvc, rc
 integer                         :: i, j, k, clb(3), cub(3)

 real(esmf_kind_r8), parameter   :: rd=287.05
 real(esmf_kind_r8), parameter   :: cp=1004.6
 real(esmf_kind_r8), parameter   :: rocp=rd/cp
 real(esmf_kind_r8), parameter   :: rocp1=rocp+1
 real(esmf_kind_r8), parameter   :: rocpr=1/rocp

 real(esmf_kind_r8), pointer     :: delp_ptr(:,:,:)
 real(esmf_kind_r8), pointer     :: pptr(:,:,:)    ! adjusted 3-d p.
 real(esmf_kind_r8), pointer     :: psptr(:,:)  ! adjusted surface p.
 real(esmf_kind_r8)              :: ak, bk
 real(esmf_kind_r8), allocatable :: pi(:,:,:)
                
 print*,"COMPUTE 3-D PRESSURE FROM ADJUSTED SURFACE PRESSURE."

 idvc = 2 ! hard wire for now.
 idsl = 2 ! hard wire for now.

 print*,"- CALL FieldGet FOR 3-D PRES."
 call ESMF_FieldGet(pres_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=pptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR DELP."
 call ESMF_FieldGet(delp_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=delp_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR SURFACE PRESSURE AFTER ADJUSTMENT"
 call ESMF_FieldGet(ps_target_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 allocate(pi(clb(1):cub(1),clb(2):cub(2),1:levp_target))

 if(idvc.eq.2) then
   do k=1,levp_target
     ak = vcoord_target(k,1) 
     bk = vcoord_target(k,2) 
     do i= clb(1), cub(1)
       do j= clb(2), cub(2)
         pi(i,j,k) = ak + bk*psptr(i,j)
       enddo
     enddo
   enddo 
   do k=1,lev_target
     do i= clb(1), cub(1)
       do j= clb(2), cub(2)
         delp_ptr(i,j,k) = pi(i,j,k) - pi(i,j,k+1)
       enddo
     enddo
   enddo
 else 
   print*,'routine only works for idvc 2'
   call mpi_abort
 endif

 if(idsl.eq.2) then 
   do k=1,lev_target
     do i= clb(1), cub(1)
       do j= clb(2), cub(2)
         pptr(i,j,k) = (pi(i,j,k)+pi(i,j,k+1))/2.0
       enddo
     enddo
   enddo 
 else 
   do k=1,lev_target
     do i= clb(1), cub(1)
       do j= clb(2), cub(2)
         pptr(i,j,k) = ((pi(i,j,k)**rocp1-pi(i,j,k+1)**rocp1)/        &
                        (rocp1*(pi(i,j,k)-pi(i,j,k+1))))**rocpr        
       enddo
     enddo
   enddo 
 endif 

 deallocate(pi)

 if (localpet == 0) then
    print*,'new pres ',pptr(clb(1),clb(2),:)
    print*,'delp     ',delp_ptr(clb(1),clb(2),:)
 endif

 end subroutine newpr1 

 subroutine newps(localpet)

!$$$  subprogram documentation block
!
! subprogram:    newps       compute new surface pressure
!   prgmmr: iredell          org: w/nmc23     date: 92-10-31
!
! abstract: computes a new surface pressure given a new orography.
!   the new pressure is computed assuming a hydrostatic balance
!   and a constant temperature lapse rate.  below ground, the
!   lapse rate is assumed to be -6.5 k/km.
!
! program history log:
!   91-10-31  mark iredell
!   2018-apr  adapt for fv3. george gayno
!
!c$$$

 implicit none

 integer, intent(in) :: localpet
 integer                         :: i, j, k
 integer                         :: clb(3), cub(3), ls, rc

 real(esmf_kind_r8), pointer     :: pptr(:,:,:)
 real(esmf_kind_r8), pointer     :: psptr(:,:)
 real(esmf_kind_r8), pointer     :: psnewptr(:,:)  ! adjusted surface p.
 real(esmf_kind_r8), pointer     :: tptr(:,:,:)
 real(esmf_kind_r8), pointer     :: qptr(:,:,:)
 real(esmf_kind_r8), pointer     :: zsptr(:,:)
 real(esmf_kind_r8), pointer     :: zsnewptr(:,:)
 real(esmf_kind_r8), allocatable :: zu(:,:)
 real(esmf_kind_r8), parameter   :: beta=-6.5E-3
 real(esmf_kind_r8), parameter   :: epsilon=1.E-9
 real(esmf_kind_r8), parameter   :: g=9.80665
 real(esmf_kind_r8), parameter   :: rd=287.05
 real(esmf_kind_r8), parameter   :: rv=461.50
 real(esmf_kind_r8), parameter   :: gor=g/rd
 real(esmf_kind_r8), parameter   :: fv=rv/rd-1.
 real(esmf_kind_r8)              :: ftv, fgam, apu, fz0
 real(esmf_kind_r8)              :: atvu, atv, fz1, fp0
 real(esmf_kind_r8)              :: apd, azd, agam, azu
 real(esmf_kind_r8)              :: atvd, fp1, gamma, pu
 real(esmf_kind_r8)              :: tvu, pd, tvd
 real(esmf_kind_r8)              :: at, aq, ap, az

 ftv(at,aq)=at*(1+fv*aq)
 fgam(apu,atvu,apd,atvd)=-gor*log(atvd/atvu)/log(apd/apu)
 fz0(ap,atv,azd,apd)=azd+atv/gor*log(apd/ap)
 fz1(ap,atv,azd,apd,agam)=azd-atv/agam*((apd/ap)**(-agam/gor)-1)
 fp0(az,azu,apu,atvu)=apu*exp(-gor/atvu*(az-azu))
 fp1(az,azu,apu,atvu,agam)=apu*(1+agam/atvu*(az-azu))**(-gor/agam)

 print*,"- ADJUST SURFACE PRESSURE FOR NEW TERRAIN."

 print*,"- CALL FieldGet FOR 3-D PRES."
 call ESMF_FieldGet(pres_b4adj_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=pptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 if(localpet==0) then
   print*,'old pres ',pptr(clb(1),clb(2),:)
 endif

 print*,"- CALL FieldGet FOR TEMPERATURE"
 call ESMF_FieldGet(temp_b4adj_target_grid, &
                    farrayPtr=tptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR SPECIFIC HUMIDITY"
 call ESMF_FieldGet(spec_hum_b4adj_target_grid, &
                    farrayPtr=qptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR SURFACE PRESSURE BEFORE ADJUSTMENT"
 call ESMF_FieldGet(ps_b4adj_target_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort


 print*,"- CALL FieldGet FOR SURFACE PRESSURE AFTER ADJUSTMENT"
 call ESMF_FieldGet(ps_target_grid, &
                    farrayPtr=psnewptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR OLD TERRAIN"
 call ESMF_FieldGet(terrain_interp_to_target_grid, &
                    farrayPtr=zsptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call mpi_abort

 print*,"- CALL FieldGet FOR NEW TERRAIN"
 call ESMF_FieldGet(terrain_target_grid, &
                    farrayPtr=zsnewptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 allocate(zu(clb(1):cub(1),clb(2):cub(2)))

!-----------------------------------------------------------------------------------
! Note, this routine was adapted from the spectral GFS which labeled the lowest
! model layer as '1'.  
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
! Compute surface pressure below the original ground.
!-----------------------------------------------------------------------------------

 ls=0
 k=1
 gamma=beta
 do i=clb(1), cub(1)
 do j=clb(2), cub(2)
   pu=pptr(i,j,k)
   tvu=ftv(tptr(i,j,k),qptr(i,j,k))
   zu(i,j)=fz1(pu,tvu,zsptr(i,j),psptr(i,j),gamma)
   if(zsnewptr(i,j).le.zu(i,j)) then
     pu=pptr(i,j,k)
     tvu=ftv(tptr(i,j,k),qptr(i,j,k))
     if(abs(gamma).gt.epsilon) then
       psnewptr(i,j)=fp1(zsnewptr(i,j),zu(i,j),pu,tvu,gamma)
     else
       psnewptr(i,j)=fp0(zsnewptr(i,j),zu(i,j),pu,tvu)
     endif
   else
     psnewptr(i,j)=0
     ls=ls+1
   endif
 enddo
 enddo

!-----------------------------------------------------------------------------------
! Compute surface pressure above the original ground.
!-----------------------------------------------------------------------------------

 do k=2,cub(3)
   if(ls.gt.0) then
     do i=clb(1),cub(1)
     do j=clb(2),cub(2)
       if(psnewptr(i,j).eq.0) then
         pu=pptr(i,j,k)
         tvu=ftv(tptr(i,j,k),qptr(i,j,k))
         pd=pptr(i,j,k-1)
         tvd=ftv(tptr(i,j,k-1),qptr(i,j,k-1))
         gamma=fgam(pu,tvu,pd,tvd)
         if(abs(gamma).gt.epsilon) then
           zu(i,j)=fz1(pu,tvu,zu(i,j),pd,gamma)
         else
           zu(i,j)=fz0(pu,tvu,zu(i,j),pd)
         endif
         if(zsnewptr(i,j).le.zu(i,j)) then
           if(abs(gamma).gt.epsilon) then
             psnewptr(i,j)=fp1(zsnewptr(i,j),zu(i,j),pu,tvu,gamma)
           else
             psnewptr(i,j)=fp0(zsnewptr(i,j),zu(i,j),pu,tvu)
           endif
           ls=ls-1
         endif
       endif
     enddo
     enddo
   endif
 enddo

!-----------------------------------------------------------------------------------
! Compute surface pressure over the top.
!-----------------------------------------------------------------------------------

 if(ls.gt.0) then
   k=cub(3)
   gamma=0
   do i=clb(1),cub(1)
   do j=clb(2),cub(2)
     if(psnewptr(i,j).eq.0) then
       pu=pptr(i,j,k)
       tvu=ftv(tptr(i,j,k),qptr(i,j,k))
       psnewptr(i,j)=fp0(zsnewptr(i,j),zu(i,j),pu,tvu)
     endif
   enddo
   enddo
 endif

 deallocate(zu)

 if (localpet == 0) then
!  do i=clb(1),cub(1)
!  do j=clb(2),cub(2)
   do i=clb(1),clb(1)
   do j=clb(2),clb(2)
     print*,'sfcp adjust ',(zsnewptr(i,j)-zsptr(i,j)), psptr(i,j),psnewptr(i,j)
   enddo
   enddo
 endif

 end subroutine newps

 subroutine read_vcoord_info

!---------------------------------------------------------------------------------
! Read vertical coordinate information.
!---------------------------------------------------------------------------------

 implicit none

 integer                    :: istat, n, k

 print*
 print*,"OPEN VERTICAL COORD FILE: ", trim(vcoord_file_target_grid)
 open(14, file=trim(vcoord_file_target_grid), form='formatted', iostat=istat)
 if (istat /= 0) then
   print*,"FATAL ERROR OPENING FILE. ISTAT IS: ", istat
   call mpi_abort
 endif

 read(14, *, iostat=istat) nvcoord_target, lev_target
 if (istat /= 0) then
   print*,"FATAL ERROR READING FILE HEADER. ISTAT IS: ",istat
   call mpi_abort
 endif

 levp_target = lev_target + 1

 allocate(vcoord_target(levp_target, nvcoord_target))
 read(14, *, iostat=istat) ((vcoord_target(n,k), k=1,nvcoord_target), n=1,levp_target)
 if (istat /= 0) then
   print*,"FATAL ERROR READING FILE. ISTAT IS: ",istat
   call errexit(6)
 endif

 print*
 do k = 1, levp_target
   print*,'VCOORD FOR LEV ', k, 'IS: ', vcoord_target(k,:)
 enddo

 close(14)

 end subroutine read_vcoord_info

 subroutine vintg(localpet)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:    VINTG       VERTICALLY INTERPOLATE UPPER-AIR FIELDS
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
!
! ABSTRACT: VERTICALLY INTERPOLATE UPPER-AIR FIELDS.
!   WIND, TEMPERATURE, HUMIDITY AND OTHER TRACERS ARE INTERPOLATED.
!   THE INTERPOLATION IS CUBIC LAGRANGIAN IN LOG PRESSURE
!   WITH A MONOTONIC CONSTRAINT IN THE CENTER OF THE DOMAIN.
!   IN THE OUTER INTERVALS IT IS LINEAR IN LOG PRESSURE.
!   OUTSIDE THE DOMAIN, FIELDS ARE GENERALLY HELD CONSTANT,
!   EXCEPT FOR TEMPERATURE AND HUMIDITY BELOW THE INPUT DOMAIN,
!   WHERE THE TEMPERATURE LAPSE RATE IS HELD FIXED AT -6.5 K/KM AND
!   THE RELATIVE HUMIDITY IS HELD CONSTANT.
!
! PROGRAM HISTORY LOG:
!   91-10-31  MARK IREDELL
!
! USAGE:    CALL VINTG(IM,KM1,KM2,NT,P1,U1,V1,T1,Q1,P2,
!    &                 U2,V2,T2,Q2)
!   INPUT ARGUMENT LIST:
!     IM           INTEGER NUMBER OF POINTS TO COMPUTE
!     KM1          INTEGER NUMBER OF INPUT LEVELS
!     KM2          INTEGER NUMBER OF OUTPUT LEVELS
!     NT           INTEGER NUMBER OF TRACERS
!     P1           REAL (IM,KM1) INPUT PRESSURES
!                  ORDERED FROM BOTTOM TO TOP OF ATMOSPHERE
!     U1           REAL (IM,KM1) INPUT ZONAL WIND
!     V1           REAL (IM,KM1) INPUT MERIDIONAL WIND
!     T1           REAL (IM,KM1) INPUT TEMPERATURE (K)
!     Q1           REAL (IM,KM1,NT) INPUT TRACERS (HUMIDITY FIRST)
!     P2           REAL (IM,KM2) OUTPUT PRESSURES
!   OUTPUT ARGUMENT LIST:
!     U2           REAL (IM,KM2) OUTPUT ZONAL WIND
!     V2           REAL (IM,KM2) OUTPUT MERIDIONAL WIND
!     T2           REAL (IM,KM2) OUTPUT TEMPERATURE (K)
!     Q2           REAL (IM,KM2,NT) OUTPUT TRACERS (HUMIDITY FIRST)
!
! SUBPROGRAMS CALLED:
!   TERP3        CUBICALLY INTERPOLATE IN ONE DIMENSION
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN
!
!C$$$
      IMPLICIT NONE

      include 'mpif.h'

       integer, intent(in) :: localpet

      INTEGER      :: IM, KM1, KM2, NT

      REAL, PARAMETER          :: DLTDZ=-6.5E-3*287.05/9.80665
      REAL, PARAMETER          :: DLPVDRT=-2.5E6/461.50

      INTEGER                  :: I, J, K, CLB(3), CUB(3), RC

      REAL                     :: DZ
        
 real(esmf_kind_r8), pointer     :: p1ptr(:,:,:)    ! input pressure
 real(esmf_kind_r8), pointer     :: p2ptr(:,:,:)    ! output pressure
 real(esmf_kind_r8), pointer     :: dzdt1ptr(:,:,:)    ! input vvel
 real(esmf_kind_r8), pointer     :: dzdt2ptr(:,:,:)    ! output vvel
 real(esmf_kind_r8), pointer     :: t1ptr(:,:,:)    ! input temperature
 real(esmf_kind_r8), pointer     :: t2ptr(:,:,:)    ! output temperature
 real(esmf_kind_r8), pointer     :: q1ptr(:,:,:)    ! input spec hum
 real(esmf_kind_r8), pointer     :: q2ptr(:,:,:)    ! output spec hum
 real(esmf_kind_r8), pointer     :: wind1ptr(:,:,:,:)  ! input wind (x,y,z components)
 real(esmf_kind_r8), pointer     :: wind2ptr(:,:,:,:)  ! input wind (x,y,z components)
 real(esmf_kind_r8), pointer     :: clw1ptr(:,:,:)  ! input wind clw
 real(esmf_kind_r8), pointer     :: clw2ptr(:,:,:)  ! output wind clw
 real(esmf_kind_r8), pointer     :: o31ptr(:,:,:)  ! input ozone
 real(esmf_kind_r8), pointer     :: o32ptr(:,:,:)  ! output ozone
 real(esmf_kind_r8), pointer     :: rwmr1ptr(:,:,:)  ! input rain water
 real(esmf_kind_r8), pointer     :: rwmr2ptr(:,:,:)  ! output rain water
 real(esmf_kind_r8), pointer     :: snmr1ptr(:,:,:)  ! input snow water
 real(esmf_kind_r8), pointer     :: snmr2ptr(:,:,:)  ! output snow water
 real(esmf_kind_r8), pointer     :: icmr1ptr(:,:,:)  ! input ice water
 real(esmf_kind_r8), pointer     :: icmr2ptr(:,:,:)  ! output ice water
 real(esmf_kind_r8), pointer     :: grle1ptr(:,:,:)  ! input graupel
 real(esmf_kind_r8), pointer     :: grle2ptr(:,:,:)  ! output graupel
 real(esmf_kind_r8), pointer     :: clda1ptr(:,:,:)  ! input cld_amt
 real(esmf_kind_r8), pointer     :: clda2ptr(:,:,:)  ! output cloud_amt

      REAL,ALLOCATABLE         :: Z1(:,:,:), Z2(:,:,:)
      REAL,ALLOCATABLE         :: C1(:,:,:,:),C2(:,:,:,:)
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE LOG PRESSURE INTERPOLATING COORDINATE
!  AND COPY INPUT WIND, TEMPERATURE, HUMIDITY AND OTHER TRACERS
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 print*,"- VERTICALY INTERPOLATE FIELDS."

 print*,"- CALL FieldGet FOR 3-D PRES."
 call ESMF_FieldGet(pres_b4adj_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=p1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR 3-D ADJUSTED PRESS"
 call ESMF_FieldGet(pres_target_grid, &
                    farrayPtr=p2ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR 3-D TEMP."
 call ESMF_FieldGet(temp_b4adj_target_grid, &
                    farrayPtr=t1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR 3-D ADJUSTED TEMP."
 call ESMF_FieldGet(temp_target_grid, &
                    farrayPtr=t2ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR 3-D Q."
 call ESMF_FieldGet(spec_hum_b4adj_target_grid, &
                    farrayPtr=q1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR 3-D ADJUSTED Q."
 call ESMF_FieldGet(spec_hum_target_grid, &
                    farrayPtr=q2ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR 3-D WIND."
 call ESMF_FieldGet(wind_b4adj_target_grid, &
                    farrayPtr=wind1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR 3-D ADJUSTED WIND."
 call ESMF_FieldGet(wind_target_grid, &
                    farrayPtr=wind2ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR VERTICAL VELOCITY."
 call ESMF_FieldGet(dzdt_b4adj_target_grid, &
                    farrayPtr=dzdt1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR ADJUSTED VERTICAL VELOCITY."
 call ESMF_FieldGet(dzdt_target_grid, &
                    farrayPtr=dzdt2ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR CLOUD LIQ WATER."
 call ESMF_FieldGet(liq_wat_b4adj_target_grid, &
                    farrayPtr=clw1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR ADJUSTED CLOUD LIQ WATER."
 call ESMF_FieldGet(liq_wat_target_grid, &
                    farrayPtr=clw2ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR OZONE."
 call ESMF_FieldGet(o3mr_b4adj_target_grid, &
                    farrayPtr=o31ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR AJUSTED OZONE."
 call ESMF_FieldGet(o3mr_target_grid, &
                    farrayPtr=o32ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort






 if (gfdl_mp) then

   print*,"- CALL FieldGet FOR RAIN WATER."
   call ESMF_FieldGet(rwmr_b4adj_target_grid, &
                      farrayPtr=rwmr1ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
           call mpi_abort

   print*,"- CALL FieldGet FOR ADJUSTED RAIN WATER."
   call ESMF_FieldGet(rwmr_target_grid, &
                      farrayPtr=rwmr2ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

   print*,"- CALL FieldGet FOR SNOW WATER."
   call ESMF_FieldGet(snmr_b4adj_target_grid, &
                      farrayPtr=snmr1ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
           call mpi_abort

   print*,"- CALL FieldGet FOR ADJUSTED SNOW WATER."
   call ESMF_FieldGet(snmr_target_grid, &
                      farrayPtr=snmr2ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

   print*,"- CALL FieldGet FOR ICE WATER."
   call ESMF_FieldGet(icmr_b4adj_target_grid, &
                      farrayPtr=icmr1ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
           call mpi_abort

   print*,"- CALL FieldGet FOR ADJUSTED ICE WATER."
   call ESMF_FieldGet(icmr_target_grid, &
                      farrayPtr=icmr2ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

   print*,"- CALL FieldGet FOR GRAUPEL."
   call ESMF_FieldGet(grle_b4adj_target_grid, &
                      farrayPtr=grle1ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
           call mpi_abort

   print*,"- CALL FieldGet FOR ADJUSTED GRAUPEL."
   call ESMF_FieldGet(grle_target_grid, &
                      farrayPtr=grle2ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

   print*,"- CALL FieldGet FOR Cloud amt."
   call ESMF_FieldGet(clda_b4adj_target_grid, &
                      farrayPtr=clda1ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
           call mpi_abort

   print*,"- CALL FieldGet FOR ADJUSTED CLOUD."
   call ESMF_FieldGet(clda_target_grid, &
                      farrayPtr=clda2ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort


 endif

 ALLOCATE(Z1(CLB(1):CUB(1),CLB(2):CUB(2),LEV_INPUT))
 ALLOCATE(Z2(CLB(1):CUB(1),CLB(2):CUB(2),LEV_TARGET))

! Note the 'z' component of the horizontal wind will be treated as a
! tracer.  So add one extra third dimension to these 3-d arrays.

 ALLOCATE(C1(CLB(1):CUB(1),CLB(2):CUB(2),LEV_INPUT,NTRACER_INPUT+5))
 ALLOCATE(C2(CLB(1):CUB(1),CLB(2):CUB(2),LEV_TARGET,NTRACER_TARGET+5))

 DO K=1,LEV_INPUT
   DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)
     Z1(I,J,K)   = -LOG(P1PTR(I,J,K))
     C1(I,J,K,1) =  WIND1PTR(I,J,K,1)  ! 'x' component of horiz wind
     C1(I,J,K,2) =  WIND1PTR(I,J,K,2)  ! 'y' component of horiz wind
     C1(I,J,K,3) =  DZDT1PTR(I,J,K)
     C1(I,J,K,4) =  T1PTR(I,J,K)
     C1(I,J,K,5) =  Q1PTR(I,J,K)
   ENDDO
   ENDDO
!  if(localpet==0) print*,'before vintg p/t:',k,p1ptr(1,1,k),c1(1,1,k,4)
 ENDDO

! hardwire tracers to three for now (plus one extra for 'z' wind).

 DO K=1,LEV_INPUT
   DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)
     C1(I,J,K,6) = CLW1PTR(I,J,K)
     C1(I,J,K,7) = O31PTR(I,J,K)
     C1(I,J,K,8) = WIND1PTR(I,J,K,3)  ! 'z' component of horiz wind
   ENDDO
   ENDDO
 ENDDO

 IF (GFDL_MP) THEN
   DO K=1,LEV_INPUT
     DO I=CLB(1),CUB(1)
     DO J=CLB(2),CUB(2)
       C1(I,J,K,9) =  RWMR1PTR(I,J,K)
       C1(I,J,K,10) = SNMR1PTR(I,J,K)
       C1(I,J,K,11) = ICMR1PTR(I,J,K)
       C1(I,J,K,12) = GRLE1PTR(I,J,K)
       C1(I,J,K,13) = CLDA1PTR(I,J,K)
     ENDDO
     ENDDO
   ENDDO
 ENDIF

 DO K=1,LEV_TARGET
   DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)
     Z2(I,J,K) = -LOG(P2PTR(I,J,K))
   ENDDO
   ENDDO
 ENDDO

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  PERFORM LAGRANGIAN ONE-DIMENSIONAL INTERPOLATION
!  THAT IS 4TH-ORDER IN INTERIOR, 2ND-ORDER IN OUTSIDE INTERVALS
!  AND 1ST-ORDER FOR EXTRAPOLATION.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 IM = (CUB(1)-CLB(1)+1) * (CUB(2)-CLB(2)+1)
 KM1= LEV_INPUT
 KM2= LEV_TARGET
 NT=  NTRACER_TARGET + 1 ! treat 'z' wind as tracer.

 CALL TERP3(IM,1,1,1,1,4+NT,(IM*KM1),(IM*KM2), &
            KM1,IM,IM,Z1,C1,KM2,IM,IM,Z2,C2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COPY OUTPUT WIND, TEMPERATURE, HUMIDITY AND OTHER TRACERS
!  EXCEPT BELOW THE INPUT DOMAIN, LET TEMPERATURE INCREASE WITH A FIXED
!  LAPSE RATE AND LET THE RELATIVE HUMIDITY REMAIN CONSTANT.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 DO K=1,LEV_TARGET
   DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)
     WIND2PTR(I,J,K,1)=C2(I,J,K,1)
     WIND2PTR(I,J,K,2)=C2(I,J,K,2)
     DZDT2PTR(I,J,K)=C2(I,J,K,3)
     DZ=Z2(I,J,K)-Z1(I,J,1)
     IF(DZ.GE.0) THEN
       T2PTR(I,J,K)=C2(I,J,K,4)
       Q2PTR(I,J,K)=C2(I,J,K,5)
     ELSE
       T2PTR(I,J,K)=T1PTR(I,J,1)*EXP(DLTDZ*DZ)
       Q2PTR(I,J,K)=Q1PTR(I,J,1)*EXP(DLPVDRT*(1/T2PTR(I,J,K)-1/T1PTR(I,J,1))-DZ)
    ENDIF
   ENDDO
   ENDDO
!  if (localpet==0) print*,'after  vintg p/t:',k,p2ptr(1,1,k),t2ptr(1,1,k)
 ENDDO

 DO K=1,LEV_TARGET
   DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)
     CLW2PTR(I,J,K)=C2(I,J,K,6)
     O32PTR(I,J,K)=C2(I,J,K,7)
     WIND2PTR(I,J,K,3)=C2(I,J,K,8)
   ENDDO
   ENDDO
 ENDDO

 IF (GFDL_MP) THEN
   DO K=1,LEV_TARGET
     DO I=CLB(1),CUB(1)
     DO J=CLB(2),CUB(2)
       RWMR2PTR(I,J,K) = C2(I,J,K,9)
       SNMR2PTR(I,J,K) = C2(I,J,K,10)
       ICMR2PTR(I,J,K) = C2(I,J,K,11)
       GRLE2PTR(I,J,K) = C2(I,J,K,12)
       CLDA2PTR(I,J,K) = C2(I,J,K,13)
     ENDDO
     ENDDO
   ENDDO
 ENDIF

 DEALLOCATE (Z1, Z2, C1, C2)

 end subroutine vintg

 subroutine terp3(IM,IXZ1,IXQ1,IXZ2,IXQ2,NM,NXQ1,NXQ2,             &
                  KM1,KXZ1,KXQ1,Z1,Q1,KM2,KXZ2,KXQ2,Z2,Q2)      
!$$$  SUBPROGRAM DOCUMENTATION BLOCK                                    
!                                                                       
! SUBPROGRAM:    TERP3       CUBICALLY INTERPOLATE IN ONE DIMENSION     
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 98-05-01            
!                                                                       
! ABSTRACT: INTERPOLATE FIELD(S) IN ONE DIMENSION ALONG THE COLUMN(S).  
!   THE INTERPOLATION IS CUBIC LAGRANGIAN WITH A MONOTONIC CONSTRAINT   
!   IN THE CENTER OF THE DOMAIN.  IN THE OUTER INTERVALS IT IS LINEAR.  
!   OUTSIDE THE DOMAIN, FIELDS ARE HELD CONSTANT.                       
!                                                                       
! PROGRAM HISTORY LOG:                                                  
!   98-05-01  MARK IREDELL                                              
! 1999-01-04  IREDELL  USE ESSL SEARCH                                  
!                                                                       
! USAGE:    CALL TERP3(IM,IXZ1,IXQ1,IXZ2,IXQ2,NM,NXQ1,NXQ2,             
!    &                 KM1,KXZ1,KXQ1,Z1,Q1,KM2,KXZ2,KXQ2,Z2,Q2,J2)      
!   INPUT ARGUMENT LIST:                                                
!     IM           INTEGER NUMBER OF COLUMNS                            
!     IXZ1         INTEGER COLUMN SKIP NUMBER FOR Z1                    
!     IXQ1         INTEGER COLUMN SKIP NUMBER FOR Q1                    
!     IXZ2         INTEGER COLUMN SKIP NUMBER FOR Z2                    
!     IXQ2         INTEGER COLUMN SKIP NUMBER FOR Q2                    
!     NM           INTEGER NUMBER OF FIELDS PER COLUMN                  
!     NXQ1         INTEGER FIELD SKIP NUMBER FOR Q1                     
!     NXQ2         INTEGER FIELD SKIP NUMBER FOR Q2                     
!     KM1          INTEGER NUMBER OF INPUT POINTS                       
!     KXZ1         INTEGER POINT SKIP NUMBER FOR Z1                     
!     KXQ1         INTEGER POINT SKIP NUMBER FOR Q1                     
!     Z1           REAL (1+(IM-1)*IXZ1+(KM1-1)*KXZ1)                    
!                  INPUT COORDINATE VALUES IN WHICH TO INTERPOLATE      
!                  (Z1 MUST BE STRICTLY MONOTONIC IN EITHER DIRECTION)  
!     Q1           REAL (1+(IM-1)*IXQ1+(KM1-1)*KXQ1+(NM-1)*NXQ1)        
!                  INPUT FIELDS TO INTERPOLATE                          
!     KM2          INTEGER NUMBER OF OUTPUT POINTS                      
!     KXZ2         INTEGER POINT SKIP NUMBER FOR Z2                     
!     KXQ2         INTEGER POINT SKIP NUMBER FOR Q2                     
!     Z2           REAL (1+(IM-1)*IXZ2+(KM2-1)*KXZ2)                    
!                  OUTPUT COORDINATE VALUES TO WHICH TO INTERPOLATE     
!                  (Z2 NEED NOT BE MONOTONIC)                           
!                                                                       
!   OUTPUT ARGUMENT LIST:                                               
!     Q2           REAL (1+(IM-1)*IXQ2+(KM2-1)*KXQ2+(NM-1)*NXQ2)        
!                  OUTPUT INTERPOLATED FIELDS                           
!     J2           REAL (1+(IM-1)*IXQ2+(KM2-1)*KXQ2+(NM-1)*NXQ2)        
!                  OUTPUT INTERPOLATED FIELDS CHANGE WRT Z2             
!                                                                       
! SUBPROGRAMS CALLED:                                                   
!   RSEARCH      SEARCH FOR A SURROUNDING REAL INTERVAL                 
!                                                                       
! ATTRIBUTES:                                                           
!   LANGUAGE: FORTRAN                                                   
!                                                                       
!C$$$                                                                   
      IMPLICIT NONE 
      INTEGER IM,IXZ1,IXQ1,IXZ2,IXQ2,NM,NXQ1,NXQ2 
      INTEGER KM1,KXZ1,KXQ1,KM2,KXZ2,KXQ2 
      INTEGER I,K1,K2,N 
      REAL Z1(1+(IM-1)*IXZ1+(KM1-1)*KXZ1) 
      REAL Q1(1+(IM-1)*IXQ1+(KM1-1)*KXQ1+(NM-1)*NXQ1) 
      REAL Z2(1+(IM-1)*IXZ2+(KM2-1)*KXZ2) 
      REAL Q2(1+(IM-1)*IXQ2+(KM2-1)*KXQ2+(NM-1)*NXQ2) 
!     REAL J2(1+(IM-1)*IXQ2+(KM2-1)*KXQ2+(NM-1)*NXQ2) 
      REAL FFA(IM),FFB(IM),FFC(IM),FFD(IM) 
      REAL GGA(IM),GGB(IM),GGC(IM),GGD(IM) 
      INTEGER K1S(IM,KM2) 
      REAL Z1A,Z1B,Z1C,Z1D,Q1A,Q1B,Q1C,Q1D,Z2S,Q2S
!     REAL J2S 


      print*,'terp3 im   ',im
      print*,'terp3 ixz1 ', ixz1
      print*,'terp3 ixq1 ', ixq1
      print*,'terp3 ixz2 ', ixz2
      print*,'terp3 ixq2 ', ixq2
      print*,'terp3 nm   ', nm
      print*,'terp3 nxq1 ', nxq1
      print*,'terp3 nxq2 ', nxq2
      print*,'terp3 km1  ', km1
      print*,'terp3 kxz1 ', kxz1
      print*,'terp3 kxq1 ', kxq1
      print*,'terp3 z1   ', z1(1), sizeof(z1)
      print*,'terp3 q1   ', q1(1), sizeof(q1)
      print*,'terp3 km2  ', km2
      print*,'terp3 kxz2 ', kxz2
      print*,'terp3 kxq2 ', kxq2
      print*,'terp3 z2   ', z2(1), sizeof(z2)
      print*,'terp3 q2   ', q2(1), sizeof(q2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!  FIND THE SURROUNDING INPUT INTERVAL FOR EACH OUTPUT POINT.           
      CALL RSEARCH(IM,KM1,IXZ1,KXZ1,Z1,KM2,IXZ2,KXZ2,Z2,1,IM,K1S) 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!  GENERALLY INTERPOLATE CUBICALLY WITH MONOTONIC CONSTRAINT            
!  FROM TWO NEAREST INPUT POINTS ON EITHER SIDE OF THE OUTPUT POINT,    
!  BUT WITHIN THE TWO EDGE INTERVALS INTERPOLATE LINEARLY.              
!  KEEP THE OUTPUT FIELDS CONSTANT OUTSIDE THE INPUT DOMAIN.            
                                                                        
!!!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(IM,IXZ1,IXQ1,IXZ2), &
!!!!$OMP& SHARED(IXQ2,NM,NXQ1,NXQ2,KM1,KXZ1,KXQ1,Z1,Q1,KM2,KXZ2), &
!!!!$OMP& SHARED(KXQ2,Z2,Q2,J2,K1S)
                                                                        
      DO K2=1,KM2 
        DO I=1,IM 
          K1=K1S(I,K2) 
          IF(K1.EQ.1.OR.K1.EQ.KM1-1) THEN 
            Z2S=Z2(1+(I-1)*IXZ2+(K2-1)*KXZ2) 
            Z1A=Z1(1+(I-1)*IXZ1+(K1-1)*KXZ1) 
            Z1B=Z1(1+(I-1)*IXZ1+(K1+0)*KXZ1) 
            FFA(I)=(Z2S-Z1B)/(Z1A-Z1B) 
            FFB(I)=(Z2S-Z1A)/(Z1B-Z1A) 
            GGA(I)=1/(Z1A-Z1B) 
            GGB(I)=1/(Z1B-Z1A) 
          ELSEIF(K1.GT.1.AND.K1.LT.KM1-1) THEN 
            Z2S=Z2(1+(I-1)*IXZ2+(K2-1)*KXZ2) 
            Z1A=Z1(1+(I-1)*IXZ1+(K1-2)*KXZ1) 
            Z1B=Z1(1+(I-1)*IXZ1+(K1-1)*KXZ1) 
            Z1C=Z1(1+(I-1)*IXZ1+(K1+0)*KXZ1) 
            Z1D=Z1(1+(I-1)*IXZ1+(K1+1)*KXZ1) 
            FFA(I)=(Z2S-Z1B)/(Z1A-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1A-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1A-Z1D)                                  
            FFB(I)=(Z2S-Z1A)/(Z1B-Z1A)*                                 &
                   (Z2S-Z1C)/(Z1B-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1B-Z1D)                                  
            FFC(I)=(Z2S-Z1A)/(Z1C-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1C-Z1B)*                                 &
                   (Z2S-Z1D)/(Z1C-Z1D)                                  
            FFD(I)=(Z2S-Z1A)/(Z1D-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1D-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1D-Z1C)                                  
            GGA(I)=        1/(Z1A-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1A-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1A-Z1D)+                                 &
                   (Z2S-Z1B)/(Z1A-Z1B)*                                 &
                           1/(Z1A-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1A-Z1D)+                                 &
                   (Z2S-Z1B)/(Z1A-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1A-Z1C)*                                 &
                           1/(Z1A-Z1D)                                  
            GGB(I)=        1/(Z1B-Z1A)*                                 &
                   (Z2S-Z1C)/(Z1B-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1B-Z1D)+                                 &
                   (Z2S-Z1A)/(Z1B-Z1A)*                                 &
                           1/(Z1B-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1B-Z1D)+                                 &
                   (Z2S-Z1A)/(Z1B-Z1A)*                                 &
                   (Z2S-Z1C)/(Z1B-Z1C)*                                 &
                           1/(Z1B-Z1D)                                  
            GGC(I)=        1/(Z1C-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1C-Z1B)*                                 &
                   (Z2S-Z1D)/(Z1C-Z1D)+                                 &
                   (Z2S-Z1A)/(Z1C-Z1A)*                                 &
                           1/(Z1C-Z1B)*                                 &
                   (Z2S-Z1D)/(Z1C-Z1D)+                                 &
                   (Z2S-Z1A)/(Z1C-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1C-Z1B)*                                 &
                           1/(Z1C-Z1D)                                  
            GGD(I)=        1/(Z1D-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1D-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1D-Z1C)+                                 &
                   (Z2S-Z1A)/(Z1D-Z1A)*                                 &
                           1/(Z1D-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1D-Z1C)+                                 &
                   (Z2S-Z1A)/(Z1D-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1D-Z1B)*                                 &
                           1/(Z1D-Z1C)                                  
          ENDIF 
        ENDDO 
!  INTERPOLATE.                                                         
        DO N=1,NM 
          DO I=1,IM 
            K1=K1S(I,K2) 
            IF(K1.EQ.0) THEN 
              Q2S=Q1(1+(I-1)*IXQ1+(N-1)*NXQ1) 
!             J2S=0 
            ELSEIF(K1.EQ.KM1) THEN 
              Q2S=Q1(1+(I-1)*IXQ1+(KM1-1)*KXQ1+(N-1)*NXQ1) 
!             J2S=0 
            ELSEIF(K1.EQ.1.OR.K1.EQ.KM1-1) THEN 
              Q1A=Q1(1+(I-1)*IXQ1+(K1-1)*KXQ1+(N-1)*NXQ1) 
              Q1B=Q1(1+(I-1)*IXQ1+(K1+0)*KXQ1+(N-1)*NXQ1) 
              Q2S=FFA(I)*Q1A+FFB(I)*Q1B 
!             J2S=GGA(I)*Q1A+GGB(I)*Q1B 
            ELSE 
              Q1A=Q1(1+(I-1)*IXQ1+(K1-2)*KXQ1+(N-1)*NXQ1) 
              Q1B=Q1(1+(I-1)*IXQ1+(K1-1)*KXQ1+(N-1)*NXQ1) 
              Q1C=Q1(1+(I-1)*IXQ1+(K1+0)*KXQ1+(N-1)*NXQ1) 
              Q1D=Q1(1+(I-1)*IXQ1+(K1+1)*KXQ1+(N-1)*NXQ1) 
              Q2S=FFA(I)*Q1A+FFB(I)*Q1B+FFC(I)*Q1C+FFD(I)*Q1D 
!             J2S=GGA(I)*Q1A+GGB(I)*Q1B+GGC(I)*Q1C+GGD(I)*Q1D 
              IF(Q2S.LT.MIN(Q1B,Q1C)) THEN 
                Q2S=MIN(Q1B,Q1C) 
!               J2S=0 
              ELSEIF(Q2S.GT.MAX(Q1B,Q1C)) THEN 
                Q2S=MAX(Q1B,Q1C) 
!               J2S=0 
              ENDIF 
            ENDIF 
            Q2(1+(I-1)*IXQ2+(K2-1)*KXQ2+(N-1)*NXQ2)=Q2S 
!           J2(1+(I-1)*IXQ2+(K2-1)*KXQ2+(N-1)*NXQ2)=J2S 
          ENDDO 
        ENDDO 
      ENDDO 
!!!!!$OMP END PARALLEL DO                                                   

 end subroutine terp3 

 subroutine rsearch(IM,KM1,IXZ1,KXZ1,Z1,KM2,IXZ2,KXZ2,Z2,IXL2,KXL2,&
                         L2)                                            
!$$$  SUBPROGRAM DOCUMENTATION BLOCK                                    
!                                                                       
! SUBPROGRAM:    RSEARCH     SEARCH FOR A SURROUNDING REAL INTERVAL     
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 98-05-01            
!                                                                       
! ABSTRACT: THIS SUBPROGRAM SEARCHES MONOTONIC SEQUENCES OF REAL NUMBERS
!   FOR INTERVALS THAT SURROUND A GIVEN SEARCH SET OF REAL NUMBERS.     
!   THE SEQUENCES MAY BE MONOTONIC IN EITHER DIRECTION; THE REAL NUMBERS
!   MAY BE SINGLE OR DOUBLE PRECISION; THE INPUT SEQUENCES AND SETS     
!   AND THE OUTPUT LOCATIONS MAY BE ARBITRARILY DIMENSIONED.            
!                                                                       
! PROGRAM HISTORY LOG:                                                  
! 1999-01-05  MARK IREDELL                                              
!                                                                       
! USAGE:    CALL RSEARCH(IM,KM1,IXZ1,KXZ1,Z1,KM2,IXZ2,KXZ2,Z2,IXL2,KXL2,
!    &                   L2)                                            
!   INPUT ARGUMENT LIST:                                                
!     IM           INTEGER NUMBER OF SEQUENCES TO SEARCH                
!     KM1          INTEGER NUMBER OF POINTS IN EACH SEQUENCE            
!     IXZ1         INTEGER SEQUENCE SKIP NUMBER FOR Z1                  
!     KXZ1         INTEGER POINT SKIP NUMBER FOR Z1                     
!     Z1           REAL (1+(IM-1)*IXZ1+(KM1-1)*KXZ1)                    
!                  SEQUENCE VALUES TO SEARCH                            
!                  (Z1 MUST BE MONOTONIC IN EITHER DIRECTION)           
!     KM2          INTEGER NUMBER OF POINTS TO SEARCH FOR               
!                  IN EACH RESPECTIVE SEQUENCE                          
!     IXZ2         INTEGER SEQUENCE SKIP NUMBER FOR Z2                  
!     KXZ2         INTEGER POINT SKIP NUMBER FOR Z2                     
!     Z2           REAL (1+(IM-1)*IXZ2+(KM2-1)*KXZ2)                    
!                  SET OF VALUES TO SEARCH FOR                          
!                  (Z2 NEED NOT BE MONOTONIC)                           
!     IXL2         INTEGER SEQUENCE SKIP NUMBER FOR L2                  
!     KXL2         INTEGER POINT SKIP NUMBER FOR L2                     
!                                                                       
!   OUTPUT ARGUMENT LIST:                                               
!     L2           INTEGER (1+(IM-1)*IXL2+(KM2-1)*KXL2)                 
!                  INTERVAL LOCATIONS HAVING VALUES FROM 0 TO KM1       
!                  (Z2 WILL BE BETWEEN Z1(L2) AND Z1(L2+1))             
!                                                                       
! SUBPROGRAMS CALLED:                                                   
!   SBSRCH       ESSL BINARY SEARCH                                     
!   DBSRCH       ESSL BINARY SEARCH                                     
!                                                                       
! REMARKS:                                                              
!   IF THE ARRAY Z1 IS DIMENSIONED (IM,KM1), THEN THE SKIP NUMBERS ARE  
!   IXZ1=1 AND KXZ1=IM; IF IT IS DIMENSIONED (KM1,IM), THEN THE SKIP    
!   NUMBERS ARE IXZ1=KM1 AND KXZ1=1; IF IT IS DIMENSIONED (IM,JM,KM1),  
!   THEN THE SKIP NUMBERS ARE IXZ1=1 AND KXZ1=IM*JM; ETCETERA.          
!   SIMILAR EXAMPLES APPLY TO THE SKIP NUMBERS FOR Z2 AND L2.           
!                                                                       
!   RETURNED VALUES OF 0 OR KM1 INDICATE THAT THE GIVEN SEARCH VALUE    
!   IS OUTSIDE THE RANGE OF THE SEQUENCE.                               
!                                                                       
!   IF A SEARCH VALUE IS IDENTICAL TO ONE OF THE SEQUENCE VALUES        
!   THEN THE LOCATION RETURNED POINTS TO THE IDENTICAL VALUE.           
!   IF THE SEQUENCE IS NOT STRICTLY MONOTONIC AND A SEARCH VALUE IS     
!   IDENTICAL TO MORE THAN ONE OF THE SEQUENCE VALUES, THEN THE         
!   LOCATION RETURNED MAY POINT TO ANY OF THE IDENTICAL VALUES.         
!                                                                       
!   TO BE EXACT, FOR EACH I FROM 1 TO IM AND FOR EACH K FROM 1 TO KM2,  
!   Z=Z2(1+(I-1)*IXZ2+(K-1)*KXZ2) IS THE SEARCH VALUE AND               
!   L=L2(1+(I-1)*IXL2+(K-1)*KXL2) IS THE LOCATION RETURNED.             
!   IF L=0, THEN Z IS LESS THAN THE START POINT Z1(1+(I-1)*IXZ1)        
!   FOR ASCENDING SEQUENCES (OR GREATER THAN FOR DESCENDING SEQUENCES). 
!   IF L=KM1, THEN Z IS GREATER THAN OR EQUAL TO THE END POINT          
!   Z1(1+(I-1)*IXZ1+(KM1-1)*KXZ1) FOR ASCENDING SEQUENCES               
!   (OR LESS THAN OR EQUAL TO FOR DESCENDING SEQUENCES).                
!   OTHERWISE Z IS BETWEEN THE VALUES Z1(1+(I-1)*IXZ1+(L-1)*KXZ1) AND   
!   Z1(1+(I-1)*IXZ1+(L-0)*KXZ1) AND MAY EQUAL THE FORMER.               
!                                                                       
! ATTRIBUTES:                                                           
!   LANGUAGE: FORTRAN                                                   
!                                                                       
!C$$$                                                                   
      IMPLICIT NONE 
      INTEGER,INTENT(IN):: IM,KM1,IXZ1,KXZ1,KM2,IXZ2,KXZ2,IXL2,KXL2 
      REAL,INTENT(IN):: Z1(1+(IM-1)*IXZ1+(KM1-1)*KXZ1) 
      REAL,INTENT(IN):: Z2(1+(IM-1)*IXZ2+(KM2-1)*KXZ2) 
      INTEGER,INTENT(OUT):: L2(1+(IM-1)*IXL2+(KM2-1)*KXL2) 
      INTEGER I,K2,L 
      REAL Z 
!C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!C  FIND THE SURROUNDING INPUT INTERVAL FOR EACH OUTPUT POINT.          
      DO I=1,IM 
        IF(Z1(1+(I-1)*IXZ1).LE.Z1(1+(I-1)*IXZ1+(KM1-1)*KXZ1)) THEN 
!C  INPUT COORDINATE IS MONOTONICALLY ASCENDING.                        
          DO K2=1,KM2 
            Z=Z2(1+(I-1)*IXZ2+(K2-1)*KXZ2) 
            L=0 
            DO 
              IF(Z.LT.Z1(1+(I-1)*IXZ1+L*KXZ1)) EXIT 
              L=L+1 
              IF(L.EQ.KM1) EXIT 
            ENDDO 
            L2(1+(I-1)*IXL2+(K2-1)*KXL2)=L 
          ENDDO 
        ELSE 
!C  INPUT COORDINATE IS MONOTONICALLY DESCENDING.                       
          DO K2=1,KM2 
            Z=Z2(1+(I-1)*IXZ2+(K2-1)*KXZ2) 
            L=0 
            DO 
              IF(Z.GT.Z1(1+(I-1)*IXZ1+L*KXZ1)) EXIT 
              L=L+1 
              IF(L.EQ.KM1) EXIT 
            ENDDO 
            L2(1+(I-1)*IXL2+(K2-1)*KXL2)=L 
          ENDDO 
        ENDIF 
      ENDDO 
                                                                        
 end subroutine rsearch 

 subroutine compute_zh

 implicit none 

 integer                         :: i,j,k, rc, clb(2), cub(2)

 real, allocatable               :: pe0(:), pn0(:)
 real(esmf_kind_r8), pointer     :: psptr(:,:)
 real(esmf_kind_r8), pointer     :: zhsfcptr(:,:)
 real(esmf_kind_r8), pointer     :: zhptr(:,:,:)
 real(esmf_kind_r8), pointer     :: tptr(:,:,:)
 real(esmf_kind_r8), pointer     :: qptr(:,:,:)
 real                            :: ak, bk, zvir, grd
 real, parameter                 :: grav  = 9.80665 
 real, parameter                 :: rdgas = 287.05 
 real, parameter                 :: rvgas = 461.50 

 print*,"- COMPUTE HEIGHT"

 print*,"- CALL FieldGet FOR SURFACE PRESSURE"
 call ESMF_FieldGet(ps_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR TERRAIN HEIGHT"
 call ESMF_FieldGet(terrain_target_grid, &
                    farrayPtr=zhsfcptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR HEIGHT"
 call ESMF_FieldGet(zh_target_grid, &
                    farrayPtr=zhptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR TEMPERATURE"
 call ESMF_FieldGet(temp_target_grid, &
                    farrayPtr=tptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 print*,"- CALL FieldGet FOR SPECIFIC HUMIDITY"
 call ESMF_FieldGet(spec_hum_target_grid, &
                    farrayPtr=qptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call mpi_abort

 grd = grav/rdgas 
 zvir = rvgas/rdgas - 1. 
                                                                        
 allocate(pe0(levp_target))
 allocate(pn0(levp_target))

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)

   do k = 1, levp_target
     ak = vcoord_target(k,1)
     ak = max(ak, 1.e-9)
     bk = vcoord_target(k,2)

     pe0(k) = ak + bk*psptr(i,j)
     pn0(k) = log(pe0(k))
   enddo

   zhptr(i,j,1) = zhsfcptr(i,j)

   do k = 2, levp_target
     zhptr(i,j,k) = zhptr(i,j,k-1)+tptr(i,j,k-1)*(1.+zvir*qptr(i,j,k-1))*     &
              (pn0(k-1)-pn0(k))/grd
   enddo

 enddo
 enddo

 deallocate(pe0, pn0)

 end subroutine compute_zh 

 subroutine cleanup_target_atm_b4adj_data

 implicit none

 integer                     :: rc

 print*,"- DESTROY TARGET GRID ATMOSPHERIC BEFORE ADJUSTMENT FIELDS."

 call ESMF_FieldDestroy(wind_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(dzdt_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(liq_wat_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(o3mr_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(ps_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(pres_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(spec_hum_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(temp_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(terrain_interp_to_target_grid, rc=rc)
 if (gfdl_mp) then
   call ESMF_FieldDestroy(clda_b4adj_target_grid, rc=rc)
   call ESMF_FieldDestroy(grle_b4adj_target_grid, rc=rc)
   call ESMF_FieldDestroy(icmr_b4adj_target_grid, rc=rc)
   call ESMF_FieldDestroy(snmr_b4adj_target_grid, rc=rc)
   call ESMF_FieldDestroy(rwmr_b4adj_target_grid, rc=rc)
 endif

 end subroutine cleanup_target_atm_b4adj_data

 subroutine cleanup_target_atm_data

 implicit none

 integer                     :: rc

 print*,"- DESTROY TARGET GRID ATMOSPHERIC FIELDS."

 call ESMF_FieldDestroy(delp_target_grid, rc=rc)
 call ESMF_FieldDestroy(dzdt_target_grid, rc=rc)
 call ESMF_FieldDestroy(liq_wat_target_grid, rc=rc)
 call ESMF_FieldDestroy(o3mr_target_grid, rc=rc)
 call ESMF_FieldDestroy(ps_target_grid, rc=rc)
 call ESMF_FieldDestroy(pres_target_grid, rc=rc)
 call ESMF_FieldDestroy(spec_hum_target_grid, rc=rc)
 call ESMF_FieldDestroy(temp_target_grid, rc=rc)
 call ESMF_FieldDestroy(u_s_target_grid, rc=rc)
 call ESMF_FieldDestroy(v_s_target_grid, rc=rc)
 call ESMF_FieldDestroy(wind_target_grid, rc=rc)
 call ESMF_FieldDestroy(wind_s_target_grid, rc=rc)
 call ESMF_FieldDestroy(wind_w_target_grid, rc=rc)
 call ESMF_FieldDestroy(u_w_target_grid, rc=rc)
 call ESMF_FieldDestroy(v_w_target_grid, rc=rc)
 call ESMF_FieldDestroy(zh_target_grid, rc=rc)
 if (gfdl_mp) then
   call ESMF_FieldDestroy(clda_target_grid, rc=rc)
   call ESMF_FieldDestroy(grle_target_grid, rc=rc)
   call ESMF_FieldDestroy(icmr_target_grid, rc=rc)
   call ESMF_FieldDestroy(snmr_target_grid, rc=rc)
   call ESMF_FieldDestroy(rwmr_target_grid, rc=rc)
 endif

 end subroutine cleanup_target_atm_data

 end module atmosphere
