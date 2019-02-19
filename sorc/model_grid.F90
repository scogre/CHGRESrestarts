 module model_grid

 use esmf

 implicit none

 private

 character(len=5), allocatable, public  :: tiles_target_grid(:)

 integer, parameter, public             :: lsoil_target = 4 ! # soil layers
 integer, public                        :: i_input, j_input
 integer, public                        :: ip1_input, jp1_input
 integer, public                        :: i_target, j_target
 integer, public                        :: ip1_target, jp1_target
 integer, public                        :: num_tiles_input_grid
 integer, public                        :: num_tiles_target_grid

 type(esmf_grid),  public               :: input_grid
 type(esmf_grid),  public               :: target_grid

 type(esmf_field),  public              :: latitude_input_grid
 type(esmf_field),  public              :: longitude_input_grid
 type(esmf_field),  public              :: latitude_s_input_grid
 type(esmf_field),  public              :: longitude_s_input_grid
 type(esmf_field),  public              :: latitude_w_input_grid
 type(esmf_field),  public              :: longitude_w_input_grid
 type(esmf_field),  public              :: latitude_c_input_grid
 type(esmf_field),  public              :: longitude_c_input_grid
 type(esmf_field),  public              :: terrain_input_grid

 type(esmf_field),  public              :: landmask_target_grid
 type(esmf_field),  public              :: latitude_target_grid
 type(esmf_field),  public              :: latitude_s_target_grid
 type(esmf_field),  public              :: latitude_w_target_grid
 type(esmf_field),  public              :: latitude_c_target_grid

 type(esmf_field),  public              :: longitude_target_grid
 type(esmf_field),  public              :: longitude_s_target_grid
 type(esmf_field),  public              :: longitude_w_target_grid
 type(esmf_field),  public              :: longitude_c_target_grid

 type(esmf_field),  public              :: seamask_target_grid
 type(esmf_field),  public              :: terrain_target_grid

 public :: define_target_grid
 public :: define_input_grid
 public :: cleanup_input_target_grid_data

 contains

 subroutine define_input_grid(localpet, npets)

 use netcdf
 use program_setup, only       : mosaic_file_input_grid,  &
                                 orog_dir_input_grid, &
                                 orog_files_input_grid

 implicit none

 include 'mpif.h'

 character(len=500)           :: the_file

 integer, intent(in)          :: localpet, npets

 integer                      :: id_tiles, id_dim, tile
 integer                      :: extra, error, ncid
 integer, allocatable         :: decomptile(:,:)

 integer(esmf_kind_i8), allocatable    :: landmask_one_tile(:,:)

 real(esmf_kind_r8), allocatable       :: latitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_s_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_w_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_c_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_s_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_w_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_c_one_tile(:,:)

 real(esmf_kind_r8), allocatable       :: terrain_one_tile(:,:)

 print*,'- OPEN INPUT GRID MOSAIC FILE: ',trim(mosaic_file_input_grid)
 error=nf90_open(trim(mosaic_file_input_grid),nf90_nowrite,ncid)
 if (error /= nf90_noerr) goto 910

 print*,"- READ NUMBER OF TILES"
 error=nf90_inq_dimid(ncid, 'ntiles', id_tiles)
 if (error /= nf90_noerr) goto 910
 error=nf90_inquire_dimension(ncid,id_tiles,len=num_tiles_input_grid)
 if (error /= nf90_noerr) goto 910

 error = nf90_close(ncid)

 print*,'- NUMBER OF TILES, INPUT MODEL GRID IS ', num_tiles_input_grid

 if (mod(npets,num_tiles_input_grid) /= 0) then
   print*,'- FATAL ERROR: MUST RUN THIS PROGRAM WITH A TASK COUNT THAT'
   print*,'- IS A MULTIPLE OF THE NUMBER OF TILES.'
   goto 910
 endif

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 extra = npets / num_tiles_input_grid

 allocate(decomptile(2,num_tiles_input_grid))

 do tile = 1, num_tiles_input_grid
   decomptile(:,tile)=(/1,extra/)
 enddo

 print*,"- CALL GridCreateMosaic FOR INPUT MODEL GRID"
 input_grid = ESMF_GridCreateMosaic(filename=trim(mosaic_file_input_grid), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, &
                                                   ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(orog_dir_input_grid), &
                                  rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

!-----------------------------------------------------------------------
! Read the mask, terrain and lat/lons.
!-----------------------------------------------------------------------

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE."
 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE."
 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE_S."
 latitude_s_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="input_grid_latitude_s", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE_S."
 longitude_s_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="input_grid_longitude_s", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE_W."
 latitude_w_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="input_grid_latitude_w", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE_W."
 longitude_w_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="input_grid_longitude_w", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE_C."
 latitude_c_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CORNER, &
                                   name="input_grid_latitude_c", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE_C."
 longitude_c_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CORNER, &
                                   name="input_grid_longitude_c", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910



 print*,"- CALL FieldCreate FOR INPUT GRID TERRAIN."
 terrain_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_terrain", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 the_file = trim(orog_dir_input_grid) // trim(orog_files_input_grid(1))

 print*,'- OPEN FIRST INPUT GRID OROGRAPHY FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 if (error /= nf90_noerr) goto 910
 print*,"- READ GRID DIMENSIONS"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 if (error /= nf90_noerr) goto 910
 error=nf90_inquire_dimension(ncid,id_dim,len=i_input)
 if (error /= nf90_noerr) goto 910
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 if (error /= nf90_noerr) goto 910
 error=nf90_inquire_dimension(ncid,id_dim,len=j_input)
 if (error /= nf90_noerr) goto 910
 error = nf90_close(ncid)

 print*,"- I/J DIMENSIONS OF THE INPUT GRID TILES ", i_input, j_input

 ip1_input = i_input + 1
 jp1_input = j_input + 1

 if (localpet == 0) then
   allocate(terrain_one_tile(i_input,j_input))
   allocate(longitude_one_tile(i_input,j_input))
   allocate(longitude_s_one_tile(i_input,jp1_input))
   allocate(longitude_w_one_tile(ip1_input,j_input))
   allocate(longitude_c_one_tile(ip1_input,jp1_input))
   allocate(latitude_one_tile(i_input,j_input))
   allocate(latitude_s_one_tile(i_input,jp1_input))
   allocate(latitude_w_one_tile(ip1_input,j_input))
   allocate(latitude_c_one_tile(ip1_input,jp1_input))
   allocate(landmask_one_tile(i_input,j_input))
 else
   allocate(terrain_one_tile(0,0))
   allocate(longitude_one_tile(0,0))
   allocate(longitude_s_one_tile(0,0))
   allocate(longitude_w_one_tile(0,0))
   allocate(longitude_c_one_tile(0,0))
   allocate(latitude_one_tile(0,0))
   allocate(latitude_s_one_tile(0,0))
   allocate(latitude_w_one_tile(0,0))
   allocate(latitude_c_one_tile(0,0))
   allocate(landmask_one_tile(0,0))
 endif

 do tile = 1, num_tiles_input_grid
   if (localpet == 0) then
     the_file = trim(orog_dir_input_grid) // trim(orog_files_input_grid(tile))
     call get_model_mask_terrain(trim(the_file), i_input, j_input, landmask_one_tile, &
                                 terrain_one_tile)
     call get_model_latlons(mosaic_file_input_grid, orog_dir_input_grid, num_tiles_input_grid, tile, &
                            i_input, j_input, ip1_input, jp1_input, latitude_one_tile, &
                            latitude_s_one_tile, latitude_w_one_tile,latitude_c_one_tile, &
                            longitude_one_tile, longitude_s_one_tile, longitude_w_one_tile, &
                            longitude_c_one_tile)
   endif
   print*,"- CALL FieldScatter FOR INPUT GRID TERRAIN. TILE IS: ", tile
   call ESMF_FieldScatter(terrain_input_grid, terrain_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_input_grid, latitude_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_input_grid, longitude_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE_S. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_s_input_grid, latitude_s_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE_S. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_s_input_grid, longitude_s_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE_W. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_w_input_grid, latitude_w_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE_W. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_w_input_grid, longitude_w_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE_C. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_c_input_grid, latitude_c_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE_C. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_c_input_grid, longitude_c_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
 enddo

 deallocate(terrain_one_tile)
 deallocate(longitude_one_tile)
 deallocate(longitude_s_one_tile)
 deallocate(longitude_w_one_tile)
 deallocate(longitude_c_one_tile)
 deallocate(latitude_one_tile)
 deallocate(latitude_s_one_tile)
 deallocate(latitude_w_one_tile)
 deallocate(latitude_c_one_tile)
 deallocate(landmask_one_tile)

 return

 910 print*,"- FATAL ERROR. STOP"
 call mpi_abort(mpi_comm_world, 18, error)

 end subroutine define_input_grid

 subroutine define_target_grid(localpet, npets)

 use netcdf
 use program_setup, only       : mosaic_file_target_grid, &
                                 orog_dir_target_grid,    &
                                 orog_files_target_grid

 implicit none

 include 'mpif.h'

 integer, intent(in)                   :: localpet, npets

 character(len=500)                    :: the_file

 integer                               :: error, ncid, extra
 integer                               :: id_tiles
 integer                               :: id_dim, id_grid_tiles
 integer                               :: tile
 integer, allocatable                  :: decomptile(:,:)
 integer(esmf_kind_i8), allocatable    :: landmask_one_tile(:,:)
 integer(esmf_kind_i8), allocatable    :: seamask_one_tile(:,:)
 
 real(esmf_kind_r8), allocatable       :: latitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_s_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_w_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_c_one_tile(:,:)

 real(esmf_kind_r8), allocatable       :: longitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_s_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_w_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_c_one_tile(:,:)

 real(esmf_kind_r8), allocatable       :: terrain_one_tile(:,:)

 print*,'- OPEN TARGET GRID MOSAIC FILE: ',trim(mosaic_file_target_grid)
 error=nf90_open(trim(mosaic_file_target_grid),nf90_nowrite,ncid)
 if (error /= nf90_noerr) goto 910

 print*,"- READ NUMBER OF TILES"
 error=nf90_inq_dimid(ncid, 'ntiles', id_tiles)
 if (error /= nf90_noerr) goto 910
 error=nf90_inquire_dimension(ncid,id_tiles,len=num_tiles_target_grid)
 if (error /= nf90_noerr) goto 910
 error=nf90_inq_varid(ncid, 'gridtiles', id_grid_tiles)
 if (error /= nf90_noerr) goto 910
 allocate(tiles_target_grid(num_tiles_target_grid))
 tiles_target_grid="NULL"
 print*,"- READ TILE NAMES"
 error=nf90_get_var(ncid, id_grid_tiles, tiles_target_grid)
 if (error /= nf90_noerr) goto 910

 error = nf90_close(ncid)

 print*,'- NUMBER OF TILES, TARGET MODEL GRID IS ', num_tiles_target_grid

 if (mod(npets,num_tiles_target_grid) /= 0) then
   print*,'- FATAL ERROR: MUST RUN THIS PROGRAM WITH A TASK COUNT THAT'
   print*,'- IS A MULTIPLE OF THE NUMBER OF TILES.'
   goto 910
 endif

!-----------------------------------------------------------------------
! Get the model grid specs and land mask from the orography files.
!-----------------------------------------------------------------------

 the_file = trim(orog_dir_target_grid) // trim(orog_files_target_grid(1))

 print*,'- OPEN FIRST TARGET GRID OROGRAPHY FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 if (error /= nf90_noerr) goto 910
 print*,"- READ GRID DIMENSIONS"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 if (error /= nf90_noerr) goto 910
 error=nf90_inquire_dimension(ncid,id_dim,len=i_target)
 if (error /= nf90_noerr) goto 910
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 if (error /= nf90_noerr) goto 910
 error=nf90_inquire_dimension(ncid,id_dim,len=j_target)
 if (error /= nf90_noerr) goto 910
 error = nf90_close(ncid)

 print*,"- I/J DIMENSIONS OF THE TARGET GRID TILES ", i_target, j_target

 ip1_target = i_target + 1
 jp1_target = j_target + 1

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 extra = npets / num_tiles_target_grid

 allocate(decomptile(2,num_tiles_target_grid))

 do tile = 1, num_tiles_target_grid
   decomptile(:,tile)=(/1,extra/)
 enddo

 print*,"- CALL GridCreateMosaic FOR TARGET GRID"
 target_grid = ESMF_GridCreateMosaic(filename=trim(mosaic_file_target_grid), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, &
                                                   ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(orog_dir_target_grid), &
                                  rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

!-----------------------------------------------------------------------
! Set target model landmask (1 - land, 0 - not land) and 
! seamask (1 - non-land, 0 -land).  Read lat/lon on target grid.
!-----------------------------------------------------------------------

 print*,"- CALL FieldCreate FOR TARGET GRID LANDMASK."
 landmask_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_I8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_landmask", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR TARGET GRID SEAMASK."
 seamask_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_I8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_seamask", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE_S."
 latitude_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_latitude_s", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE_W."
 latitude_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_latitude_w", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
! SG SG
 print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE_C."
 latitude_c_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CORNER, &
                                   name="target_grid_latitude_c", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
! SG SG

 print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE_S."
 longitude_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_longitude_s", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE_W."
 longitude_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_longitude_w", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

! SG SG
 print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE_C."
 longitude_c_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CORNER, &
                                   name="target_grid_longitude_c", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
! SG SG

 print*,"- CALL FieldCreate FOR TARGET GRID TERRAIN."
 terrain_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_terrain", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910

 if (localpet == 0) then
   allocate(landmask_one_tile(i_target,j_target))
   allocate(seamask_one_tile(i_target,j_target))
   allocate(latitude_one_tile(i_target,j_target))
   allocate(latitude_s_one_tile(i_target,jp1_target))
   allocate(latitude_w_one_tile(ip1_target,j_target))
   allocate(latitude_c_one_tile(ip1_target,jp1_target))
   allocate(longitude_one_tile(i_target,j_target))
   allocate(longitude_s_one_tile(i_target,jp1_target))
   allocate(longitude_w_one_tile(ip1_target,j_target))
   allocate(longitude_c_one_tile(ip1_target,jp1_target))
   allocate(terrain_one_tile(i_target,j_target))
 else
   allocate(landmask_one_tile(0,0))
   allocate(seamask_one_tile(0,0))
   allocate(longitude_one_tile(0,0))
   allocate(longitude_s_one_tile(0,0))
   allocate(longitude_w_one_tile(0,0))
   allocate(longitude_c_one_tile(0,0))
   allocate(latitude_one_tile(0,0))
   allocate(latitude_s_one_tile(0,0))
   allocate(latitude_w_one_tile(0,0))
   allocate(latitude_c_one_tile(0,0))
   allocate(terrain_one_tile(0,0))
 endif

 do tile = 1, num_tiles_target_grid
   if (localpet == 0) then
     the_file = trim(orog_dir_target_grid) // trim(orog_files_target_grid(tile))
     call get_model_mask_terrain(trim(the_file), i_target, j_target, landmask_one_tile, &
                                 terrain_one_tile)
     seamask_one_tile = 0
     where(landmask_one_tile == 0) seamask_one_tile = 1
     call get_model_latlons(mosaic_file_target_grid, orog_dir_target_grid,num_tiles_target_grid, tile, &
                            i_target, j_target, ip1_target, jp1_target,latitude_one_tile, &
                            latitude_s_one_tile,latitude_w_one_tile,latitude_c_one_tile, &
                            longitude_one_tile, longitude_s_one_tile,longitude_w_one_tile, &
                            longitude_c_one_tile)
   endif
   print*,"- CALL FieldScatter FOR TARGET GRID LANDMASK. TILE IS: ", tile
   call ESMF_FieldScatter(landmask_target_grid, landmask_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID SEAMASK. TILE IS: ", tile
   call ESMF_FieldScatter(seamask_target_grid, seamask_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_target_grid, longitude_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE_S. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_s_target_grid, longitude_s_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE_W. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_w_target_grid, longitude_w_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE_C. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_c_target_grid, longitude_c_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_target_grid, latitude_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE_S. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_s_target_grid, latitude_s_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE_W. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_w_target_grid, latitude_w_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE_C. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_c_target_grid, latitude_c_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))goto 910
   print*,"- CALL FieldScatter FOR TARGET GRID TERRAIN. TILE IS: ", tile
   call ESMF_FieldScatter(terrain_target_grid, terrain_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) goto 910
 enddo

 deallocate(landmask_one_tile)
 deallocate(seamask_one_tile)
 deallocate(longitude_one_tile)
 deallocate(longitude_s_one_tile)
 deallocate(longitude_w_one_tile)
 deallocate(longitude_c_one_tile)
 deallocate(latitude_one_tile)
 deallocate(latitude_s_one_tile)
 deallocate(latitude_w_one_tile)
 deallocate(latitude_c_one_tile)
 deallocate(terrain_one_tile)

 return

 910 print*,"- FATAL ERROR. STOP"
 call mpi_abort(mpi_comm_world, 14, error)

 end subroutine define_target_grid

!-----------------------------------------------------------------------
! Read model lat/lons for a single tile from the "grid" file.
!-----------------------------------------------------------------------

 subroutine get_model_latlons(mosaic_file, orog_dir, num_tiles, tile, &
                              i_tile, j_tile, ip1_tile, jp1_tile,  &
                              latitude, latitude_s, latitude_w, &
                              latitude_c, longitude, longitude_s, &
                              longitude_w, longitude_c)

 use netcdf

 implicit none

 include "mpif.h"

 character(len=*), intent(in)      :: mosaic_file, orog_dir

 integer, intent(in)               :: num_tiles, tile
 integer, intent(in)               :: i_tile, j_tile
 integer, intent(in)               :: ip1_tile, jp1_tile

 real(esmf_kind_r8), intent(out)   :: latitude(i_tile, j_tile)
 real(esmf_kind_r8), intent(out)   :: latitude_s(i_tile, jp1_tile)
 real(esmf_kind_r8), intent(out)   :: latitude_w(ip1_tile, j_tile)
 real(esmf_kind_r8), intent(out)   :: latitude_c(ip1_tile, jp1_tile)

 real(esmf_kind_r8), intent(out)   :: longitude(i_tile, j_tile)
 real(esmf_kind_r8), intent(out)   :: longitude_s(i_tile, jp1_tile)
 real(esmf_kind_r8), intent(out)   :: longitude_w(ip1_tile, j_tile)
 real(esmf_kind_r8), intent(out)   :: longitude_c(ip1_tile, jp1_tile)

 character(len=25)                 :: grid_files(num_tiles)
 character(len=255)                :: grid_file

 integer                           :: error, id_var, ncid
 integer                           :: id_dim, nxp, nyp, i, j, ii, jj

 real, allocatable                 :: tmpvar(:,:)

 print*,"- READ MODEL GRID FILE"

 print*,'- OPEN MOSAIC FILE: ', trim(mosaic_file)
 error=nf90_open(trim(mosaic_file), nf90_nowrite, ncid)
 if (error /= nf90_noerr) goto 900

 print*,"- READ GRID FILE NAMES"
 error=nf90_inq_varid(ncid, 'gridfiles', id_var)
 if (error /= nf90_noerr) goto 900
 error=nf90_get_var(ncid, id_var, grid_files)
 if (error /= nf90_noerr) goto 900

 error = nf90_close(ncid)

 grid_file = trim(orog_dir) // trim(grid_files(tile))

 print*,'- OPEN GRID FILE: ', trim(grid_file)
 error=nf90_open(trim(grid_file), nf90_nowrite, ncid)
 if (error /= nf90_noerr) goto 900

 print*,'- READ NXP ID'
 error=nf90_inq_dimid(ncid, 'nxp', id_dim)
 if (error /= nf90_noerr) goto 900

 print*,'- READ NXP'
 error=nf90_inquire_dimension(ncid,id_dim,len=nxp)
 if (error /= nf90_noerr) goto 900

 print*,'- READ NYP ID'
 error=nf90_inq_dimid(ncid, 'nyp', id_dim)
 if (error /= nf90_noerr) goto 900

 print*,'- READ NYP'
 error=nf90_inquire_dimension(ncid,id_dim,len=nyp)
 if (error /= nf90_noerr) goto 900

 if ((nxp/2 /= i_tile) .or. (nyp/2 /= j_tile)) then
   print*,'- DIMENSION MISMATCH IN GRID FILE.'
   goto 900
 endif

 allocate(tmpvar(nxp,nyp))

 print*,'- READ LONGITUDE ID'
 error=nf90_inq_varid(ncid, 'x', id_var)
 if (error /= nf90_noerr) goto 900

 print*,'- READ LONGITUDE'
 error=nf90_get_var(ncid, id_var, tmpvar)
 if (error /= nf90_noerr) goto 900

 do j = 1, j_tile
 do i = 1, i_tile
   ii = 2*i
   jj = 2*j
   longitude(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 do j = 1, jp1_tile
 do i = 1, i_tile
   ii = 2*i
   jj = (2*j) - 1
   longitude_s(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 do j = 1, j_tile
 do i = 1, ip1_tile
   ii = (2*i) - 1
   jj = 2*j
   longitude_w(i,j) = tmpvar(ii,jj)
 enddo
 enddo
!SG SG
 do j = 1, jp1_tile
 do i = 1, ip1_tile
   ii = (2*i) - 1
   jj = (2*j) - 1
   longitude_c(i,j) = tmpvar(ii,jj)
 enddo
 enddo
!SG SG



 print*,'- READ LATITUDE ID'
 error=nf90_inq_varid(ncid, 'y', id_var)
 if (error /= nf90_noerr) goto 900

 print*,'- READ LATIITUDE'
 error=nf90_get_var(ncid, id_var, tmpvar)
 if (error /= nf90_noerr) goto 900

 do j = 1, j_tile
 do i = 1, i_tile
   ii = 2*i
   jj = 2*j
   latitude(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 do j = 1, jp1_tile
 do i = 1, i_tile
   ii = 2*i
   jj = (2*j) - 1
   latitude_s(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 do j = 1, j_tile
 do i = 1, ip1_tile
   ii = (2*i) - 1
   jj = 2*j
   latitude_w(i,j) = tmpvar(ii,jj)
 enddo
 enddo

!SG SG
 do j = 1, jp1_tile
 do i = 1, ip1_tile
   ii = (2*i) - 1
   jj = (2*j) - 1
   latitude_c(i,j) = tmpvar(ii,jj)
 enddo
 enddo
!SG SG

 deallocate(tmpvar)

 error = nf90_close(ncid)

 return

 900 print*,"- FATAL ERROR. STOP"
 call mpi_abort(mpi_comm_world, 24, error)

 end subroutine get_model_latlons

!-----------------------------------------------------------------------
! Read the model land mask and terrain for a single tile.
!-----------------------------------------------------------------------

 subroutine get_model_mask_terrain(orog_file, idim, jdim, mask, terrain)

 use esmf
 use netcdf

 implicit none

 include "mpif.h"

 character(len=*), intent(in)       :: orog_file

 integer, intent(in)                :: idim, jdim
 integer(esmf_kind_i8), intent(out) :: mask(idim,jdim)

 real(esmf_kind_i8), intent(out)    :: terrain(idim,jdim)

 integer                            :: error, lat, lon
 integer                            :: ncid, id_dim, id_var

 real(kind=4), allocatable          :: dummy(:,:)

 print*,"- READ MODEL LAND MASK FILE"

 print*,'- OPEN LAND MASK FILE: ', orog_file
 error=nf90_open(orog_file,nf90_nowrite,ncid)
 if (error /= nf90_noerr) goto 900

 print*,"- READ I-DIMENSION"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 if (error /= nf90_noerr) goto 900
 error=nf90_inquire_dimension(ncid,id_dim,len=lon)
 if (error /= nf90_noerr) goto 900

 print*,"- READ J-DIMENSION"
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 if (error /= nf90_noerr) goto 900
 error=nf90_inquire_dimension(ncid,id_dim,len=lat)
 if (error /= nf90_noerr) goto 900

 print*,"- I/J DIMENSIONS: ", lon, lat

 if ((lon /= idim) .or. (lat /= jdim)) then
   print*,'- MISMATCH IN DIMENSIONS.'
   goto 900
 endif

 allocate(dummy(idim,jdim))

 print*,"- READ LAND MASK"
 error=nf90_inq_varid(ncid, 'slmsk', id_var)
 if (error /= nf90_noerr) goto 900
 error=nf90_get_var(ncid, id_var, dummy)
 if (error /= nf90_noerr) goto 900
 mask = nint(dummy)

! print*,"- READ RAW OROGRAPHY."
! error=nf90_inq_varid(ncid, 'orog_raw', id_var)
 print*,"- READ FILT OROGRAPHY."
 error=nf90_inq_varid(ncid, 'orog_filt', id_var)
 if (error /= nf90_noerr) goto 900
 error=nf90_get_var(ncid, id_var, dummy)
 if (error /= nf90_noerr) goto 900
 terrain = dummy

 error = nf90_close(ncid)

 deallocate (dummy)

 return

 900 print*,"- FATAL ERROR.  STOP."
 call mpi_abort(mpi_comm_world, 45, error)

 end subroutine get_model_mask_terrain

 subroutine cleanup_input_target_grid_data

 implicit none

 integer                                :: rc

 print*,"- DESTROY MODEL DATA."

 call ESMF_FieldDestroy(terrain_input_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(latitude_s_input_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(latitude_w_input_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(latitude_c_input_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort


 call ESMF_FieldDestroy(longitude_s_input_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(longitude_w_input_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(longitude_c_input_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort


 call ESMF_FieldDestroy(landmask_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(latitude_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(latitude_s_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(latitude_w_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(latitude_c_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(longitude_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(longitude_s_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(longitude_w_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(longitude_c_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(seamask_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_FieldDestroy(terrain_target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_GridDestroy(input_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 call ESMF_GridDestroy(target_grid, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) call mpi_abort

 end subroutine cleanup_input_target_grid_data

 end module model_grid
