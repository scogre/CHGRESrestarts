 program chgres

 use esmf

 use atmosphere, only          : atmosphere_driver

 use program_setup, only       : read_setup_namelist, &
                                 convert_atm, &
                                 convert_sfc

 use model_grid, only          : define_target_grid,  &
                                 define_input_grid, &
                                 cleanup_input_target_grid_data

 use surface, only             : surface_driver

 implicit none

 integer                      :: ierr, localpet, npets

 type(esmf_vm)                :: vm

 include 'mpif.h'

 call mpi_init(ierr)

 print*,"- INITIALIZE ESMF"
 call ESMF_Initialize(rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort(mpi_comm_world, 66, ierr)

 print*,"- CALL VMGetGlobal"
 call ESMF_VMGetGlobal(vm, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort(mpi_comm_world, 67, ierr)

 print*,"- CALL VMGet"
 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call mpi_abort(mpi_comm_world, 68, ierr)

 print*,'- NPETS IS  ',npets
 print*,'- LOCAL PET ',localpet

 call read_setup_namelist

 call define_input_grid(localpet, npets)

 call define_target_grid(localpet, npets)

 if (convert_atm) then

   call atmosphere_driver(localpet)

 end if

 if (convert_sfc) then

   call surface_driver(localpet)

 end if

 call cleanup_input_target_grid_data

 print*,"- CALL ESMF_finalize"
 call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)

 call mpi_finalize(ierr)

 print*,"- DONE."

 end program chgres
