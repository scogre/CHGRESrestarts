 subroutine netcdf_err( err, string )

 use netcdf

 implicit none
 integer, intent(in) :: err
 character(len=*), intent(in) :: string
 character(len=256) :: errmsg

 include "mpif.h"

 if( err.EQ.NF90_NOERR )return
 errmsg = NF90_STRERROR(err)
 print*,''
 print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
 print*,'STOP.'
 call mpi_abort(mpi_comm_world, 999)

 return
 end subroutine netcdf_err

 subroutine write_fv3_atm_header_netcdf(localpet)

 use netcdf

 use atmosphere, only : nvcoord_target, &
                        vcoord_target,  &
                        levp_target,    &
                        ntracer_target

 implicit none

 integer, intent(in) :: localpet

 character(len=13)   :: outfile

 integer             :: fsize=65536, inital = 0
 integer             :: header_buffer_val = 16384
 integer             :: error, ncid, dim_nvcoord
 integer             :: dim_levp, id_ntrac, id_vcoord

 real(kind=8), allocatable :: tmp(:,:)

 if (localpet /= 0) return

 print*,"- WRITE ATMOSPHERIC HEADER FILE."

 outfile="./gfs_ctrl.nc"

 error = nf90_create(outfile, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                     ncid, initialsize=inital, chunksize=fsize)
 call netcdf_err(error, 'CREATING FILE='//trim(outfile) )

 error = nf90_def_dim(ncid, 'nvcoord', nvcoord_target, dim_nvcoord)
 call netcdf_err(error, 'define dimension nvcoord for file='//trim(outfile) )

 error = nf90_def_dim(ncid, 'levsp', levp_target, dim_levp)
 call netcdf_err(error, 'define dimension levsp for file='//trim(outfile) )

 error = nf90_def_var(ncid, 'ntrac', nf90_int, id_ntrac)
 call netcdf_err(error, 'define var ntrac for file='//trim(outfile) )

 error = nf90_def_var(ncid, 'vcoord', nf90_double, (/dim_levp, dim_nvcoord/), id_vcoord)
 call netcdf_err(error, 'define var vcoord for file='//trim(outfile) )

 error = nf90_enddef(ncid, header_buffer_val,4,0,4)
 call netcdf_err(error, 'end meta define for file='//trim(outfile) )

 error = nf90_put_var( ncid, id_ntrac, ntracer_target)
 call netcdf_err(error, 'write var ntrac for file='//trim(outfile) )

 allocate(tmp(levp_target, nvcoord_target))
 tmp(1:levp_target,:) = vcoord_target(levp_target:1:-1,:)

 error = nf90_put_var( ncid, id_vcoord, tmp)
 call netcdf_err(error, 'write var vcoord for file='//trim(outfile) )

 deallocate(tmp)

 error = nf90_close(ncid)

 end subroutine write_fv3_atm_header_netcdf

 subroutine write_fv3_atm_bndy_data_netcdf(localpet)

!---------------------------------------------------------------------------
!
! Output data along the four halo boundaries.  The naming convention
! assumes point (1,1) is the lower left corner of the grid:
!
!          --------------- TOP ---------------
!          |                                 |
!          |                                 |
!     LEFT |                                 | RIGHT
!          |                                 |
!          |PT(1,1)                          |
!          ------------- BOTTOM --------------
!
!---------------------------------------------------------------------------

 use esmf
 use netcdf

 use atmosphere, only            : lev_target, levp_target, &
                                   dzdt_target_grid, &
                                   liq_wat_target_grid, &
                                   ps_target_grid, &
                                   o3mr_target_grid, &
                                   spec_hum_target_grid, &
                                   rwmr_target_grid, &
                                   icmr_target_grid, &
                                   snmr_target_grid, &
                                   grle_target_grid, &
                                   u_s_target_grid, &
                                   v_s_target_grid, &
                                   u_w_target_grid, &
                                   v_w_target_grid, &
                                   zh_target_grid

 use input_data, only            : flip

 use model_grid, only            : i_target, ip1_target, j_target, jp1_target

 use program_setup, only         : halo, halo_p1, gfdl_mp

 implicit none

 integer, intent(in)            :: localpet

 integer                        :: fsize=65536, inital = 0
 integer                        :: header_buffer_val = 16384
 integer                        :: ncid, error, tile, i
 integer                        :: dim_lon, dim_lat
 integer                        :: dim_lonp, dim_halo
 integer                        :: dim_halop, dim_latm
 integer                        :: dim_lev, dim_levp
 integer                        :: j_target2
 integer                        :: id_i_bottom, id_j_bottom
 integer                        :: id_i_top, id_j_top
 integer                        :: id_i_right, id_j_right
 integer                        :: id_i_left, id_j_left
 integer                        :: id_ps_bottom, id_ps_top
 integer                        :: id_ps_right, id_ps_left
 integer                        :: id_w_bottom, id_w_top
 integer                        :: id_w_right, id_w_left
 integer                        :: id_zh_bottom, id_zh_top
 integer                        :: id_zh_right, id_zh_left
 integer                        :: id_sphum_bottom, id_sphum_top
 integer                        :: id_sphum_right, id_sphum_left
 integer                        :: id_o3mr_bottom, id_o3mr_top
 integer                        :: id_o3mr_right, id_o3mr_left
 integer                        :: id_clwmr_bottom, id_clwmr_top
 integer                        :: id_clwmr_right, id_clwmr_left
 integer                        :: id_rwmr_bottom, id_rwmr_top
 integer                        :: id_rwmr_right, id_rwmr_left
 integer                        :: id_icmr_bottom, id_icmr_top
 integer                        :: id_icmr_right, id_icmr_left
 integer                        :: id_snmr_bottom, id_snmr_top
 integer                        :: id_snmr_right, id_snmr_left
 integer                        :: id_grle_bottom, id_grle_top
 integer                        :: id_grle_right, id_grle_left
 integer                        :: id_i_w_bottom, id_j_w_bottom
 integer                        :: id_i_w_top, id_j_w_top
 integer                        :: id_j_w_right, id_i_w_left
 integer                        :: id_j_w_left, id_i_w_right
 integer                        :: id_u_w_bottom, id_u_w_top
 integer                        :: id_u_w_right, id_u_w_left
 integer                        :: id_v_w_bottom, id_v_w_top
 integer                        :: id_v_w_right, id_v_w_left
 integer                        :: id_i_s_bottom, id_j_s_bottom
 integer                        :: id_i_s_top, id_j_s_top
 integer                        :: id_i_s_right, id_j_s_right
 integer                        :: id_i_s_left, id_j_s_left
 integer                        :: id_u_s_bottom, id_u_s_top
 integer                        :: id_u_s_right, id_u_s_left
 integer                        :: id_v_s_bottom, id_v_s_top
 integer                        :: id_v_s_right, id_v_s_left
 integer                        :: i_start_top, i_end_top
 integer                        :: j_start_top, j_end_top
 integer                        :: i_start_bottom, i_end_bottom
 integer                        :: j_start_bottom, j_end_bottom
 integer                        :: i_start_left, i_end_left
 integer                        :: j_start_left, j_end_left
 integer                        :: i_start_right, i_end_right
 integer                        :: j_start_right, j_end_right
 integer(kind=4), allocatable   :: idum(:)

 real(kind=4), allocatable        :: dum2d_top(:,:), dum2d_bottom(:,:)
 real(kind=4), allocatable        :: dum2d_left(:,:), dum2d_right(:,:)
 real(kind=4), allocatable        :: dum3d_top(:,:,:), dum3d_bottom(:,:,:)
 real(kind=4), allocatable        :: dum3d_left(:,:,:), dum3d_right(:,:,:)
 real(esmf_kind_r8), allocatable  :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable  :: data_one_tile_3d(:,:,:)

 print*,"- OUTPUT LATERAL BOUNDARY DATA."

 if (localpet == 0) then

!--- open the file
   error = nf90_create("./gfs.bndy.nc", IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                     ncid, initialsize=inital, chunksize=fsize)
   call netcdf_err(error, 'CREATING BNDY FILE' )

   error = nf90_def_dim(ncid, 'lon', i_target, dim_lon)
   call netcdf_err(error, 'defining lon dimension')

   j_target2 = j_target - (2*halo)
   error = nf90_def_dim(ncid, 'lat', j_target2, dim_lat)
   call netcdf_err(error, 'DEFINING LAT DIMENSION')

   error = nf90_def_dim(ncid, 'lonp', ip1_target, dim_lonp)
   call netcdf_err(error, 'DEFINING LONP DIMENSION')

   j_target2 = jp1_target - (2*halo_p1)
   error = nf90_def_dim(ncid, 'latm', j_target2, dim_latm)
   call netcdf_err(error, 'DEFINING LATM DIMENSION')

   error = nf90_def_dim(ncid, 'halo', halo, dim_halo)
   call netcdf_err(error, 'DEFINING HALO DIMENSION')

   error = nf90_def_dim(ncid, 'halop', halo_p1, dim_halop)
   call netcdf_err(error, 'DEFINING HALOP DIMENSION')

   error = nf90_def_dim(ncid, 'lev', lev_target, dim_lev)
   call netcdf_err(error, 'DEFINING LEV DIMENSION')

   error = nf90_def_dim(ncid, 'levp', levp_target, dim_levp)
   call netcdf_err(error, 'DEFINING LEVP DIMENSION')

   error = nf90_def_var(ncid, 'i_bottom', NF90_INT, &
                             (/dim_lon/), id_i_bottom)
   call netcdf_err(error, 'DEFINING I_BOTTOM')

   error = nf90_def_var(ncid, 'j_bottom', NF90_INT, &
                             (/dim_halo/), id_j_bottom)
   call netcdf_err(error, 'DEFINING J_BOTTOM')

   error = nf90_def_var(ncid, 'i_top', NF90_INT, &
                             (/dim_lon/), id_i_top)
   call netcdf_err(error, 'DEFINING I_TOP')

   error = nf90_def_var(ncid, 'j_top', NF90_INT, &
                             (/dim_halo/), id_j_top)
   call netcdf_err(error, 'DEFINING J_TOP')

   error = nf90_def_var(ncid, 'i_right', NF90_INT, &
                             (/dim_halo/), id_i_right)
   call netcdf_err(error, 'DEFINING I_RIGHT')

   error = nf90_def_var(ncid, 'j_right', NF90_INT, &
                             (/dim_lat/), id_j_right)
   call netcdf_err(error, 'DEFINING J_RIGHT')

   error = nf90_def_var(ncid, 'i_left', NF90_INT, &
                             (/dim_halo/), id_i_left)
   call netcdf_err(error, 'DEFINING I_LEFT')

   error = nf90_def_var(ncid, 'j_left', NF90_INT, &
                             (/dim_lat/), id_j_left)
   call netcdf_err(error, 'DEFINING J_LEFT')

   error = nf90_def_var(ncid, 'ps_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo/), id_ps_bottom)
   call netcdf_err(error, 'DEFINING PS_BOTTOM')

   error = nf90_def_var(ncid, 'ps_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo/), id_ps_top)
   call netcdf_err(error, 'DEFINING PS_TOP')

   error = nf90_def_var(ncid, 'ps_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat/), id_ps_right)
   call netcdf_err(error, 'DEFINING PS_RIGHT')

   error = nf90_def_var(ncid, 'ps_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat/), id_ps_left)
   call netcdf_err(error, 'DEFINING PS_LEFT')

   error = nf90_def_var(ncid, 'w_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_w_bottom)
   call netcdf_err(error, 'DEFINING W_BOTTOM')

   error = nf90_def_var(ncid, 'w_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_w_top)
   call netcdf_err(error, 'DEFINING W_TOP')

   error = nf90_def_var(ncid, 'w_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_w_right)
   call netcdf_err(error, 'DEFINING W_RIGHT')

   error = nf90_def_var(ncid, 'w_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_w_left)
   call netcdf_err(error, 'DEFINING W_LEFT')

   error = nf90_def_var(ncid, 'zh_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_levp/), id_zh_bottom)
   call netcdf_err(error, 'DEFINING ZH_BOTTOM')

   error = nf90_def_var(ncid, 'zh_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_levp/), id_zh_top)
   call netcdf_err(error, 'DEFINING ZH_TOP')

   error = nf90_def_var(ncid, 'zh_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_levp/), id_zh_right)
   call netcdf_err(error, 'DEFINING ZH_RIGHT')

   error = nf90_def_var(ncid, 'zh_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_levp/), id_zh_left)
   call netcdf_err(error, 'DEFINING ZH_LEFT')

   error = nf90_def_var(ncid, 'sphum_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_sphum_bottom)
   call netcdf_err(error, 'DEFINING SPHUM_BOTTOM')

   error = nf90_def_var(ncid, 'sphum_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_sphum_top)
   call netcdf_err(error, 'DEFINING SPHUM_TOP')

   error = nf90_def_var(ncid, 'sphum_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_sphum_right)
   call netcdf_err(error, 'DEFINING SPHUM_RIGHT')

   error = nf90_def_var(ncid, 'sphum_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_sphum_left)
   call netcdf_err(error, 'DEFINING SPHUM_LEFT')

   error = nf90_def_var(ncid, 'o3mr_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_o3mr_bottom)
   call netcdf_err(error, 'DEFINING O3MR_BOTTOM')

   error = nf90_def_var(ncid, 'o3mr_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_o3mr_top)
   call netcdf_err(error, 'DEFINING O3MR_TOP')

   error = nf90_def_var(ncid, 'o3mr_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_o3mr_right)
   call netcdf_err(error, 'DEFINING O3MR_RIGHT')

   error = nf90_def_var(ncid, 'o3mr_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_o3mr_left)
   call netcdf_err(error, 'DEFINING O3MR_LEFT')

   error = nf90_def_var(ncid, 'liq_wat_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_clwmr_bottom)
   call netcdf_err(error, 'DEFINING LIQ_WAT_BOTTOM')

   error = nf90_def_var(ncid, 'liq_wat_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_clwmr_top)
   call netcdf_err(error, 'DEFINING LIQ_WAT_TOP')

   error = nf90_def_var(ncid, 'liq_wat_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_clwmr_right)
   call netcdf_err(error, 'DEFINING LIQ_WAT_RIGHT')

   error = nf90_def_var(ncid, 'liq_wat_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_clwmr_left)
   call netcdf_err(error, 'DEFINING LIQ_WAT_RIGHT')

   if (gfdl_mp) then

     error = nf90_def_var(ncid, 'rwmr_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_rwmr_bottom)
     call netcdf_err(error, 'DEFINING RWMR_BOTTOM')

     error = nf90_def_var(ncid, 'rwmr_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_rwmr_top)
     call netcdf_err(error, 'DEFINING RWMR_TOP')

     error = nf90_def_var(ncid, 'rwmr_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_rwmr_right)
     call netcdf_err(error, 'DEFINING RWMR_RIGHT')

     error = nf90_def_var(ncid, 'rwmr_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_rwmr_left)
     call netcdf_err(error, 'DEFINING RWMR_RIGHT')

     error = nf90_def_var(ncid, 'icmr_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_icmr_bottom)
     call netcdf_err(error, 'DEFINING ICMR_BOTTOM')

     error = nf90_def_var(ncid, 'icmr_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_icmr_top)
     call netcdf_err(error, 'DEFINING ICMR_TOP')

     error = nf90_def_var(ncid, 'icmr_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_icmr_right)
     call netcdf_err(error, 'DEFINING ICMR_RIGHT')

     error = nf90_def_var(ncid, 'icmr_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_icmr_left)
     call netcdf_err(error, 'DEFINING ICMR_RIGHT')

     error = nf90_def_var(ncid, 'snmr_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_snmr_bottom)
     call netcdf_err(error, 'DEFINING SNMR_BOTTOM')

     error = nf90_def_var(ncid, 'snmr_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_snmr_top)
     call netcdf_err(error, 'DEFINING SNMR_TOP')

     error = nf90_def_var(ncid, 'snmr_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_snmr_right)
     call netcdf_err(error, 'DEFINING SNMR_RIGHT')

     error = nf90_def_var(ncid, 'snmr_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_snmr_left)
     call netcdf_err(error, 'DEFINING SNMR_RIGHT')

     error = nf90_def_var(ncid, 'grle_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_grle_bottom)
     call netcdf_err(error, 'DEFINING GRLE_BOTTOM')

     error = nf90_def_var(ncid, 'grle_top', NF90_FLOAT, &
                             (/dim_lon, dim_halo, dim_lev/), id_grle_top)
     call netcdf_err(error, 'DEFINING GRLE_TOP')

     error = nf90_def_var(ncid, 'grle_right', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_grle_right)
     call netcdf_err(error, 'DEFINING GRLE_RIGHT')

     error = nf90_def_var(ncid, 'grle_left', NF90_FLOAT, &
                             (/dim_halo, dim_lat, dim_lev/), id_grle_left)
     call netcdf_err(error, 'DEFINING GRLE_RIGHT')

   endif

   error = nf90_def_var(ncid, 'i_w_bottom', NF90_INT, &
                             (/dim_lonp/), id_i_w_bottom)
   call netcdf_err(error, 'DEFINING I_W_BOTTOM')

   error = nf90_def_var(ncid, 'j_w_bottom', NF90_INT, &
                             (/dim_halo/), id_j_w_bottom)
   call netcdf_err(error, 'DEFINING J_W_BOTTOM')

   error = nf90_def_var(ncid, 'i_w_top', NF90_INT, &
                             (/dim_lonp/), id_i_w_top)
   call netcdf_err(error, 'DEFINING I_W_TOP')

   error = nf90_def_var(ncid, 'j_w_top', NF90_INT, &
                             (/dim_halo/), id_j_w_top)
   call netcdf_err(error, 'DEFINING J_W_TOP')

   error = nf90_def_var(ncid, 'i_w_right', NF90_INT, &
                             (/dim_halop/), id_i_w_right)
   call netcdf_err(error, 'DEFINING I_W_RIGHT')

   error = nf90_def_var(ncid, 'j_w_right', NF90_INT, &
                             (/dim_lat/), id_j_w_right)
   call netcdf_err(error, 'DEFINING J_W_RIGHT')

   error = nf90_def_var(ncid, 'i_w_left', NF90_INT, &
                             (/dim_halop/), id_i_w_left)
   call netcdf_err(error, 'DEFINING I_W_LEFT')

   error = nf90_def_var(ncid, 'j_w_left', NF90_INT, &
                             (/dim_lat/), id_j_w_left)
   call netcdf_err(error, 'DEFINING J_W_LEFT')

   error = nf90_def_var(ncid, 'u_w_bottom', NF90_FLOAT, &
                             (/dim_lonp, dim_halo, dim_lev/), id_u_w_bottom)
   call netcdf_err(error, 'DEFINING U_W_BOTTOM')

   error = nf90_def_var(ncid, 'u_w_top', NF90_FLOAT, &
                             (/dim_lonp, dim_halo, dim_lev/), id_u_w_top)
   call netcdf_err(error, 'DEFINING U_W_TOP')

   error = nf90_def_var(ncid, 'u_w_right', NF90_FLOAT, &
                             (/dim_halop, dim_lat, dim_lev/), id_u_w_right)
   call netcdf_err(error, 'DEFINING U_W_RIGHT')

   error = nf90_def_var(ncid, 'u_w_left', NF90_FLOAT, &
                             (/dim_halop, dim_lat, dim_lev/), id_u_w_left)
   call netcdf_err(error, 'DEFINING U_W_LEFT')

   error = nf90_def_var(ncid, 'v_w_bottom', NF90_FLOAT, &
                             (/dim_lonp, dim_halo, dim_lev/), id_v_w_bottom)
   call netcdf_err(error, 'DEFINING V_W_BOTTOM')

   error = nf90_def_var(ncid, 'v_w_top', NF90_FLOAT, &
                             (/dim_lonp, dim_halo, dim_lev/), id_v_w_top)
   call netcdf_err(error, 'DEFINING V_W_TOP')

   error = nf90_def_var(ncid, 'v_w_right', NF90_FLOAT, &
                             (/dim_halop, dim_lat, dim_lev/), id_v_w_right)
   call netcdf_err(error, 'DEFINING V_W_RIGHT')

   error = nf90_def_var(ncid, 'v_w_left', NF90_FLOAT, &
                             (/dim_halop, dim_lat, dim_lev/), id_v_w_left)
   call netcdf_err(error, 'DEFINING V_W_LEFT')

   error = nf90_def_var(ncid, 'i_s_bottom', NF90_INT, &
                             (/dim_lon/), id_i_s_bottom)
   call netcdf_err(error, 'DEFINING I_S_BOTTOM')

   error = nf90_def_var(ncid, 'j_s_bottom', NF90_INT, &
                             (/dim_halop/), id_j_s_bottom)
   call netcdf_err(error, 'DEFINING J_S_BOTTOM')

   error = nf90_def_var(ncid, 'i_s_top', NF90_INT, &
                             (/dim_lon/), id_i_s_top)
   call netcdf_err(error, 'DEFINING I_S_TOP')

   error = nf90_def_var(ncid, 'j_s_top', NF90_INT, &
                             (/dim_halop/), id_j_s_top)
   call netcdf_err(error, 'DEFINING J_S_TOP')

   error = nf90_def_var(ncid, 'i_s_right', NF90_INT, &
                             (/dim_halo/), id_i_s_right)
   call netcdf_err(error, 'DEFINING I_S_RIGHT')

   error = nf90_def_var(ncid, 'j_s_right', NF90_INT, &
                             (/dim_latm/), id_j_s_right)
   call netcdf_err(error, 'DEFINING J_S_RIGHT')

   error = nf90_def_var(ncid, 'i_s_left', NF90_INT, &
                             (/dim_halo/), id_i_s_left)
   call netcdf_err(error, 'DEFINING I_S_LEFT')

   error = nf90_def_var(ncid, 'j_s_left', NF90_INT, &
                             (/dim_latm/), id_j_s_left)
   call netcdf_err(error, 'DEFINING J_S_LEFT')

   error = nf90_def_var(ncid, 'u_s_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halop, dim_lev/), id_u_s_bottom)
   call netcdf_err(error, 'DEFINING U_S_BOTTOM')

   error = nf90_def_var(ncid, 'u_s_top', NF90_FLOAT, &
                             (/dim_lon, dim_halop, dim_lev/), id_u_s_top)
   call netcdf_err(error, 'DEFINING U_S_TOP')

   error = nf90_def_var(ncid, 'u_s_right', NF90_FLOAT, &
                             (/dim_halo, dim_latm, dim_lev/), id_u_s_right)
   call netcdf_err(error, 'DEFINING U_S_RIGHT')

   error = nf90_def_var(ncid, 'u_s_left', NF90_FLOAT, &
                             (/dim_halo, dim_latm, dim_lev/), id_u_s_left)
   call netcdf_err(error, 'DEFINING U_S_LEFT')

   error = nf90_def_var(ncid, 'v_s_bottom', NF90_FLOAT, &
                             (/dim_lon, dim_halop, dim_lev/), id_v_s_bottom)
   call netcdf_err(error, 'DEFINING V_S_BOTTOM')

   error = nf90_def_var(ncid, 'v_s_top', NF90_FLOAT, &
                             (/dim_lon, dim_halop, dim_lev/), id_v_s_top)
   call netcdf_err(error, 'DEFINING V_S_TOP')

   error = nf90_def_var(ncid, 'v_s_right', NF90_FLOAT, &
                             (/dim_halo, dim_latm, dim_lev/), id_v_s_right)
   call netcdf_err(error, 'DEFINING V_S_RIGHT')

   error = nf90_def_var(ncid, 'v_s_left', NF90_FLOAT, &
                             (/dim_halo, dim_latm, dim_lev/), id_v_s_left)
   call netcdf_err(error, 'DEFINING V_S_LEFT')

   error = nf90_enddef(ncid, header_buffer_val,4,0,4)
   call netcdf_err(error, 'DEFINING END OF HEADER')

 endif

! Set up bounds

 i_start_top = 1
 i_end_top   = i_target
 j_start_top = j_target - halo + 1
 j_end_top   = j_target

 i_start_bottom = 1
 i_end_bottom   = i_target
 j_start_bottom = 1
 j_end_bottom   = halo

 i_start_left = 1
 i_end_left   = halo
 j_start_left = halo + 1
 j_end_left   = j_target - halo

 i_start_right = i_target - halo + 1
 i_end_right   = i_target
 j_start_right = halo + 1
 j_end_right   = j_target - halo

 if (localpet == 0) then
   allocate(idum(i_start_top:i_end_top))
   do i = i_start_top, i_end_top
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_top, idum)
   call netcdf_err(error, "WRITING I_TOP")
   deallocate(idum)
   allocate(idum(i_start_bottom:i_end_bottom))
   do i = i_start_bottom, i_end_bottom
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_bottom, idum)
   call netcdf_err(error, "WRITING I_BOTTOM")
   deallocate(idum)
   allocate(idum(i_start_left:i_end_left))
   do i = i_start_left, i_end_left
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_left, idum)
   call netcdf_err(error, "WRITING I_LEFT")
   deallocate(idum)
   allocate(idum(i_start_right:i_end_right))
   do i = i_start_right, i_end_right
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_right, idum)
   call netcdf_err(error, "WRITING I_RIGHT")
   deallocate(idum)
   allocate(idum(j_start_top:j_end_top))
   do i = j_start_top, j_end_top
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_top, idum)
   call netcdf_err(error, "WRITING J_TOP")
   deallocate(idum)
   allocate(idum(j_start_bottom:j_end_bottom))
   do i = j_start_bottom, j_end_bottom
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_bottom, idum)
   call netcdf_err(error, "WRITING J_BOTTOM")
   deallocate(idum)
   allocate(idum(j_start_left:j_end_left))
   do i = j_start_left, j_end_left
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_left, idum)
   call netcdf_err(error, "WRITING J_LEFT")
   deallocate(idum)
   allocate(idum(j_start_right:j_end_right))
   do i = j_start_right, j_end_right
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_right, idum)
   call netcdf_err(error, "WRITING J_RIGHT")
   deallocate(idum)
 endif

!  surface pressure

 if (localpet == 0) then
   allocate(data_one_tile(i_target,j_target))
   allocate(dum2d_top(i_target,halo))
   allocate(dum2d_bottom(i_target,halo))
   allocate(dum2d_left(halo, j_target-2*halo))
   allocate(dum2d_right(halo, j_target-2*halo))
 else
   allocate(data_one_tile(0,0))
   allocate(dum2d_top(0,0))
   allocate(dum2d_bottom(0,0))
   allocate(dum2d_left(0,0))
   allocate(dum2d_right(0,0))
 endif

 tile = 1

 print*,"- CALL FieldGather FOR TARGET GRID SURFACE PRESSURE"
 call ESMF_FieldGather(ps_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   dum2d_top(:,:) = data_one_tile(i_start_top:i_end_top, j_start_top:j_end_top)
   error = nf90_put_var( ncid, id_ps_top, dum2d_top)
   call netcdf_err(error, 'WRITING PS TOP' )
   dum2d_bottom(:,:) = data_one_tile(i_start_bottom:i_end_bottom, j_start_bottom:j_end_bottom)
   error = nf90_put_var( ncid, id_ps_bottom, dum2d_bottom)
   call netcdf_err(error, 'WRITING PS BOTTOM' )
   dum2d_left(:,:) = data_one_tile(i_start_left:i_end_left, j_start_left:j_end_left)
   error = nf90_put_var( ncid, id_ps_left, dum2d_left)
   call netcdf_err(error, 'WRITING PS LEFT' )
   dum2d_right(:,:) = data_one_tile(i_start_right:i_end_right, j_start_right:j_end_right)
   error = nf90_put_var( ncid, id_ps_right, dum2d_right)
   call netcdf_err(error, 'WRITING PS RIGHT' )
 endif

 deallocate(dum2d_top, dum2d_bottom, dum2d_left, dum2d_right, data_one_tile)

!  height

 if (localpet == 0) then
   allocate(data_one_tile_3d(i_target,j_target,levp_target))
   allocate(dum3d_top(i_target,halo,levp_target))
   allocate(dum3d_bottom(i_target,halo,levp_target))
   allocate(dum3d_left(halo, (j_target-2*halo), levp_target))
   allocate(dum3d_right(halo, (j_target-2*halo), levp_target))
 else
   allocate(data_one_tile_3d(0,0,0))
   allocate(dum3d_top(0,0,0))
   allocate(dum3d_bottom(0,0,0))
   allocate(dum3d_left(0,0,0))
   allocate(dum3d_right(0,0,0))
 endif

 print*,"- CALL FieldGather FOR TARGET GRID HEIGHT FOR TILE: ", tile
 call ESMF_FieldGather(zh_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, i_target, j_target, levp_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_zh_top, dum3d_top)
   call netcdf_err(error, 'WRITING ZH TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_zh_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING ZH BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_zh_left, dum3d_left)
   call netcdf_err(error, 'WRITING ZH LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_zh_right, dum3d_right)
   call netcdf_err(error, 'WRITING ZH RIGHT' )
 endif

 deallocate(dum3d_top, dum3d_bottom, dum3d_left, dum3d_right, data_one_tile_3d)
 
!  specific humidity

 if (localpet == 0) then
   allocate(data_one_tile_3d(i_target,j_target,lev_target))
   allocate(dum3d_top(i_target,halo,lev_target))
   allocate(dum3d_bottom(i_target,halo,lev_target))
   allocate(dum3d_left(halo, (j_target-2*halo), lev_target))
   allocate(dum3d_right(halo, (j_target-2*halo), lev_target))
 else
   allocate(data_one_tile_3d(0,0,0))
   allocate(dum3d_top(0,0,0))
   allocate(dum3d_bottom(0,0,0))
   allocate(dum3d_left(0,0,0))
   allocate(dum3d_right(0,0,0))
 endif

 print*,"- CALL FieldGather FOR TARGET GRID SPEC HUM FOR TILE: ", tile
 call ESMF_FieldGather(spec_hum_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, i_target, j_target, lev_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_sphum_top, dum3d_top)
   call netcdf_err(error, 'WRITING SPEC HUM TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_sphum_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING SPEC HUM BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_sphum_left, dum3d_left)
   call netcdf_err(error, 'WRITING SPEC HUM LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_sphum_right, dum3d_right)
   call netcdf_err(error, 'WRITING SPEC HUM RIGHT' )
 endif

! Ozone

 print*,"- CALL FieldGather FOR TARGET GRID OZONE FOR TILE: ", tile
 call ESMF_FieldGather(o3mr_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, i_target, j_target, lev_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_o3mr_top, dum3d_top)
   call netcdf_err(error, 'WRITING O3MR TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_o3mr_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING O3MR BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_o3mr_left, dum3d_left)
   call netcdf_err(error, 'WRITING O3MR LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_o3mr_right, dum3d_right)
   call netcdf_err(error, 'WRITING O3MR RIGHT' )
 endif

! Vertical velocity

 print*,"- CALL FieldGather FOR TARGET GRID W FOR TILE: ", tile
 call ESMF_FieldGather(dzdt_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, i_target, j_target, lev_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_w_top, dum3d_top)
   call netcdf_err(error, 'WRITING W TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_w_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING W BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_w_left, dum3d_left)
   call netcdf_err(error, 'WRITING W LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_w_right, dum3d_right)
   call netcdf_err(error, 'WRITING W RIGHT' )
 endif

! Cloud liq water

 print*,"- CALL FieldGather FOR TARGET GRID CLWMR FOR TILE: ", tile
 call ESMF_FieldGather(liq_wat_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, i_target, j_target, lev_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_clwmr_top, dum3d_top)
   call netcdf_err(error, 'WRITING CLWMR TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_clwmr_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING CLWMR BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_clwmr_left, dum3d_left)
   call netcdf_err(error, 'WRITING CLWMR LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_clwmr_right, dum3d_right)
   call netcdf_err(error, 'WRITING CLWMR RIGHT' )
 endif

 if (gfdl_mp) then

!  Rain water

   print*,"- CALL FieldGather FOR TARGET GRID RWMR FOR TILE: ", tile
   call ESMF_FieldGather(rwmr_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
     error = nf90_put_var( ncid, id_rwmr_top, dum3d_top)
     call netcdf_err(error, 'WRITING RWMR TOP' )
     dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
     error = nf90_put_var( ncid, id_rwmr_bottom, dum3d_bottom)
     call netcdf_err(error, 'WRITING RWMR BOTTOM' )
     dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
     error = nf90_put_var( ncid, id_rwmr_left, dum3d_left)
     call netcdf_err(error, 'WRITING RWMR LEFT' )
     dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
     error = nf90_put_var( ncid, id_rwmr_right, dum3d_right)
     call netcdf_err(error, 'WRITING RWMR RIGHT' )
   endif

!  Snow water

   print*,"- CALL FieldGather FOR TARGET GRID SNMR FOR TILE: ", tile
   call ESMF_FieldGather(snmr_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
     error = nf90_put_var( ncid, id_snmr_top, dum3d_top)
     call netcdf_err(error, 'WRITING SNMR TOP' )
     dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
     error = nf90_put_var( ncid, id_snmr_bottom, dum3d_bottom)
     call netcdf_err(error, 'WRITING SNMR BOTTOM' )
     dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
     error = nf90_put_var( ncid, id_snmr_left, dum3d_left)
     call netcdf_err(error, 'WRITING SNMR LEFT' )
     dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
     error = nf90_put_var( ncid, id_snmr_right, dum3d_right)
     call netcdf_err(error, 'WRITING SNMR RIGHT' )
   endif

!  Ice water

   print*,"- CALL FieldGather FOR TARGET GRID ICMR FOR TILE: ", tile
   call ESMF_FieldGather(icmr_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
     error = nf90_put_var( ncid, id_icmr_top, dum3d_top)
     call netcdf_err(error, 'WRITING ICMR TOP' )
     dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
     error = nf90_put_var( ncid, id_icmr_bottom, dum3d_bottom)
     call netcdf_err(error, 'WRITING ICMR BOTTOM' )
     dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
     error = nf90_put_var( ncid, id_icmr_left, dum3d_left)
     call netcdf_err(error, 'WRITING ICMR LEFT' )
     dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
     error = nf90_put_var( ncid, id_icmr_right, dum3d_right)
     call netcdf_err(error, 'WRITING ICMR RIGHT' )
   endif

!  Graupel

   print*,"- CALL FieldGather FOR TARGET GRID GRAUPEL FOR TILE: ", tile
   call ESMF_FieldGather(grle_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
     error = nf90_put_var( ncid, id_grle_top, dum3d_top)
     call netcdf_err(error, 'WRITING GRLE TOP' )
     dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
     error = nf90_put_var( ncid, id_grle_bottom, dum3d_bottom)
     call netcdf_err(error, 'WRITING GRLE BOTTOM' )
     dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
     error = nf90_put_var( ncid, id_grle_left, dum3d_left)
     call netcdf_err(error, 'WRITING GRLE LEFT' )
     dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
     error = nf90_put_var( ncid, id_grle_right, dum3d_right)
     call netcdf_err(error, 'WRITING GRLE RIGHT' )
   endif

 endif

 deallocate(dum3d_top, dum3d_bottom, dum3d_left, dum3d_right, data_one_tile_3d)

! Set up bounds for staggered 'S' winds

 i_start_top = 1
 i_end_top   = i_target
 j_start_top = jp1_target - halo_p1 + 1
 j_end_top   = jp1_target

 i_start_bottom = 1
 i_end_bottom   = i_target
 j_start_bottom = 1
 j_end_bottom   = halo_p1

 i_start_left = 1
 i_end_left   = halo
 j_start_left = halo_p1 + 1
 j_end_left   = jp1_target - halo_p1

 i_start_right = i_target - halo + 1
 i_end_right   = i_target
 j_start_right = halo_p1 + 1
 j_end_right   = jp1_target - halo_p1

 if (localpet == 0) then
   allocate(idum(i_start_top:i_end_top))
   do i = i_start_top, i_end_top
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_s_top, idum)
   call netcdf_err(error, "WRITING I_S_TOP")
   deallocate(idum)
   allocate(idum(i_start_bottom:i_end_bottom))
   do i = i_start_bottom, i_end_bottom
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_s_bottom, idum)
   call netcdf_err(error, "WRITING I_S_BOTTOM")
   deallocate(idum)
   allocate(idum(i_start_left:i_end_left))
   do i = i_start_left, i_end_left
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_s_left, idum)
   call netcdf_err(error, "WRITING I_S_LEFT")
   deallocate(idum)
   allocate(idum(i_start_right:i_end_right))
   do i = i_start_right, i_end_right
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_s_right, idum)
   call netcdf_err(error, "WRITING I_S_RIGHT")
   deallocate(idum)
   allocate(idum(j_start_top:j_end_top))
   do i = j_start_top, j_end_top
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_s_top, idum)
   call netcdf_err(error, "WRITING J_S_TOP")
   deallocate(idum)
   allocate(idum(j_start_bottom:j_end_bottom))
   do i = j_start_bottom, j_end_bottom
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_s_bottom, idum)
   call netcdf_err(error, "WRITING J_S_BOTTOM")
   deallocate(idum)
   allocate(idum(j_start_left:j_end_left))
   do i = j_start_left, j_end_left
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_s_left, idum)
   call netcdf_err(error, "WRITING J_S_LEFT")
   deallocate(idum)
   allocate(idum(j_start_right:j_end_right))
   do i = j_start_right, j_end_right
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_s_right, idum)
   call netcdf_err(error, "WRITING J_S_RIGHT")
   deallocate(idum)
 endif

! U-WINDS 'S'

 if (localpet == 0) then
   allocate(data_one_tile_3d(i_target,jp1_target,lev_target))
   allocate(dum3d_top(i_target,halo_p1,lev_target))
   allocate(dum3d_bottom(i_target,halo_p1,lev_target))
   allocate(dum3d_left(halo, (j_end_left-j_start_left+1), lev_target))
   allocate(dum3d_right(halo, (j_end_right-j_start_right+1), lev_target))
 else
   allocate(data_one_tile_3d(0,0,0))
   allocate(dum3d_top(0,0,0))
   allocate(dum3d_bottom(0,0,0))
   allocate(dum3d_left(0,0,0))
   allocate(dum3d_right(0,0,0))
 endif

 print*,"- CALL FieldGather FOR TARGET GRID U_S FOR TILE: ", tile
 call ESMF_FieldGather(u_s_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, i_target, jp1_target, lev_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_u_s_top, dum3d_top)
   call netcdf_err(error, 'WRITING U_S TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_u_s_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING U_S BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_u_s_left, dum3d_left)
   call netcdf_err(error, 'WRITING U_S LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_u_s_right, dum3d_right)
   call netcdf_err(error, 'WRITING U_S RIGHT' )
 endif

! V-WINDS 'S'

 print*,"- CALL FieldGather FOR TARGET GRID V_S FOR TILE: ", tile
 call ESMF_FieldGather(v_s_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, i_target, jp1_target, lev_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_v_s_top, dum3d_top)
   call netcdf_err(error, 'WRITING V_S TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_v_s_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING V_S BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_v_s_left, dum3d_left)
   call netcdf_err(error, 'WRITING V_S LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_v_s_right, dum3d_right)
   call netcdf_err(error, 'WRITING V_S RIGHT' )
 endif

 deallocate(dum3d_top, dum3d_bottom, dum3d_left, dum3d_right, data_one_tile_3d)

! Set up bounds for staggered 'W' winds

 i_start_top = 1
 i_end_top   = ip1_target
 j_start_top = j_target - halo + 1
 j_end_top   = j_target

 i_start_bottom = 1
 i_end_bottom   = ip1_target
 j_start_bottom = 1
 j_end_bottom   = halo

 i_start_left = 1
 i_end_left   = halo_p1
 j_start_left = halo_p1
 j_end_left   = j_target - halo

 i_start_right = ip1_target - halo_p1 + 1
 i_end_right   = ip1_target
 j_start_right = halo_p1
 j_end_right   = j_target - halo

 if (localpet == 0) then
   allocate(idum(i_start_top:i_end_top))
   do i = i_start_top, i_end_top
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_w_top, idum)
   call netcdf_err(error, "WRITING I_W_TOP")
   deallocate(idum)
   allocate(idum(i_start_bottom:i_end_bottom))
   do i = i_start_bottom, i_end_bottom
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_w_bottom, idum)
   call netcdf_err(error, "WRITING I_W_BOTTOM")
   deallocate(idum)
   allocate(idum(i_start_left:i_end_left))
   do i = i_start_left, i_end_left
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_w_left, idum)
   call netcdf_err(error, "WRITING I_W_LEFT")
   deallocate(idum)
   allocate(idum(i_start_right:i_end_right))
   do i = i_start_right, i_end_right
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_i_w_right, idum)
   call netcdf_err(error, "WRITING I_W_RIGHT")
   deallocate(idum)
   allocate(idum(j_start_top:j_end_top))
   do i = j_start_top, j_end_top
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_w_top, idum)
   call netcdf_err(error, "WRITING J_W_TOP")
   deallocate(idum)
   allocate(idum(j_start_bottom:j_end_bottom))
   do i = j_start_bottom, j_end_bottom
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_w_bottom, idum)
   call netcdf_err(error, "WRITING J_W_BOTTOM")
   deallocate(idum)
   allocate(idum(j_start_left:j_end_left))
   do i = j_start_left, j_end_left
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_w_left, idum)
   call netcdf_err(error, "WRITING J_W_LEFT")
   deallocate(idum)
   allocate(idum(j_start_right:j_end_right))
   do i = j_start_right, j_end_right
     idum(i) = i
   enddo
   error = nf90_put_var(ncid, id_j_w_right, idum)
   call netcdf_err(error, "WRITING J_W_RIGHT")
   deallocate(idum)
 endif

! U-WINDS 'W'

 if (localpet == 0) then
   allocate(data_one_tile_3d(ip1_target,j_target,lev_target))
   allocate(dum3d_top(ip1_target,halo,lev_target))
   allocate(dum3d_bottom(ip1_target,halo,lev_target))
   allocate(dum3d_left(halo_p1, (j_end_left-j_start_left+1), lev_target))
   allocate(dum3d_right(halo_p1, (j_end_right-j_start_right+1), lev_target))
 else
   allocate(data_one_tile_3d(0,0,0))
   allocate(dum3d_top(0,0,0))
   allocate(dum3d_bottom(0,0,0))
   allocate(dum3d_left(0,0,0))
   allocate(dum3d_right(0,0,0))
 endif

 print*,"- CALL FieldGather FOR TARGET GRID U_W FOR TILE: ", tile
 call ESMF_FieldGather(u_w_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, ip1_target, j_target, lev_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_u_w_top, dum3d_top)
   call netcdf_err(error, 'WRITING U_W TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_u_w_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING U_W BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_u_w_left, dum3d_left)
   call netcdf_err(error, 'WRITING U_W LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_u_w_right, dum3d_right)
   call netcdf_err(error, 'WRITING U_W RIGHT' )
 endif

! V-WINDS 'W'

 print*,"- CALL FieldGather FOR TARGET GRID V_W FOR TILE: ", tile
 call ESMF_FieldGather(v_w_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

 if (localpet == 0) then
   call flip(data_one_tile_3d, ip1_target, j_target, lev_target)
   dum3d_top(:,:,:) = data_one_tile_3d(i_start_top:i_end_top,j_start_top:j_end_top,:)
   error = nf90_put_var( ncid, id_v_w_top, dum3d_top)
   call netcdf_err(error, 'WRITING V_W TOP' )
   dum3d_bottom(:,:,:) = data_one_tile_3d(i_start_bottom:i_end_bottom,j_start_bottom:j_end_bottom,:)
   error = nf90_put_var( ncid, id_v_w_bottom, dum3d_bottom)
   call netcdf_err(error, 'WRITING V_W BOTTOM' )
   dum3d_left(:,:,:) = data_one_tile_3d(i_start_left:i_end_left,j_start_left:j_end_left,:)
   error = nf90_put_var( ncid, id_v_w_left, dum3d_left)
   call netcdf_err(error, 'WRITING V_W LEFT' )
   dum3d_right(:,:,:) = data_one_tile_3d(i_start_right:i_end_right,j_start_right:j_end_right,:)
   error = nf90_put_var( ncid, id_v_w_right, dum3d_right)
   call netcdf_err(error, 'WRITING V_W RIGHT' )
 endif

 deallocate(dum3d_top, dum3d_bottom, dum3d_left, dum3d_right, data_one_tile_3d)

 if (localpet == 0) error = nf90_close(ncid)

 end subroutine write_fv3_atm_bndy_data_netcdf

 subroutine write_fv3_atm_data_netcdf(localpet)

 use esmf
 use netcdf

 use program_setup, only           : halo, gfdl_mp

 use atmosphere, only              : lev_target, &
                                     levp_target, &
                                     ntracer_target, &
                                     liq_wat_target_grid, &
                                     o3mr_target_grid, &
                                     ps_target_grid, &
                                     zh_target_grid, &
                                     dzdt_target_grid, &
                                     spec_hum_target_grid, &
                                     rwmr_target_grid, &
                                     snmr_target_grid, &
                                     icmr_target_grid, &
                                     grle_target_grid, &
                                     temp_target_grid, &
                                     delp_target_grid, &
                                     u_s_target_grid,   &
                                     v_s_target_grid,   &
                                     u_w_target_grid,   &
                                     v_w_target_grid

 use model_grid, only              : num_tiles_target_grid, &
                                     i_target, j_target, &
                                     ip1_target, jp1_target, &
                                     terrain_target_grid, &
                                     longitude_target_grid, &
                                     longitude_s_target_grid, &
                                     longitude_w_target_grid, &
                                     longitude_c_target_grid, &
                                     latitude_s_target_grid, &
                                     latitude_w_target_grid, &
                                     latitude_c_target_grid, &
                                     latitude_target_grid

 use input_data, only              : flip

 implicit none

 integer, intent(in)              :: localpet

 character(len=128)               :: outfile

 integer                          :: error, ncid, tile
 integer                          :: fsize=65536, inital = 0
 integer                          :: header_buffer_val = 16384
 integer                          :: dim_lon, dim_lat, dim_time
 integer                          :: dim_lonp, dim_latp
 integer                          :: dim_lev, dim_levp, dim_ntracer
 integer                          :: id_sphum, id_o3mr, id_liq_wat
 integer                          :: id_lon, id_lat, id_ps
 integer                          :: id_lon_w, id_lat_w
 integer                          :: id_x, id_xp, id_y, id_yp
 integer                          :: id_lon_s, id_lat_s
 integer                          :: id_lon_c, id_lat_c
 integer                          :: id_w, id_zh, id_u_w, id_time
 integer                          :: id_z, id_v_w, id_u_s, id_v_s
 integer                          :: id_u, id_v
 integer                          :: id_t, id_delp, id_phis
 integer                          :: id_rwmr, id_icmr, id_snmr, id_grle
 integer                          :: i_start, i_end, j_start, j_end
 integer                          :: i_target_out, j_target_out
 integer                          :: ip1_target_out, jp1_target_out
 integer                          :: ip1_end, jp1_end
 integer                          :: i,j,k,nlev

 real, parameter                  :: gravity = 9.8
 real(kind=4), allocatable        :: lev_data(:), x_data(:), y_data(:)
 real(kind=4), allocatable        :: x_datap1(:), y_datap1(:)
 real(esmf_kind_r8), allocatable  :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable  :: data_one_tile_3d(:,:,:)
 real(esmf_kind_r8), allocatable  :: data_one_tile_ijp1(:,:)
 real(esmf_kind_r8), allocatable  :: data_one_tile_ip1j(:,:)
 real(esmf_kind_r8), allocatable  :: data_one_tile_ip1jp1(:,:)
 real(kind=4), allocatable        :: dum2d(:,:)
 real(kind=4), allocatable        :: dum2d_ijp1(:,:)
 real(kind=4), allocatable        :: dum2d_ip1j(:,:)
 real(kind=4), allocatable        :: dum2d_ip1jp1(:,:)
 real(kind=4), allocatable        :: lat_corner(:,:)
 real(kind=4), allocatable        :: lon_corner(:,:)
 real(kind=4), allocatable        :: dum3d(:,:,:)
 real(kind=4), allocatable        :: u_s(:,:,:)
 real(kind=4), allocatable        :: v_s(:,:,:)
 real(kind=4), allocatable        :: u_w(:,:,:)
 real(kind=4), allocatable        :: v_w(:,:,:)
 real(kind=4)                     :: times
 real,allocatable,dimension(:)    :: p1,p2,p3,ex,ey,e1
 real,allocatable,dimension(:,:,:) ::ud,vd,ua,va,udnew,vdnew
 real inner_prod
 real, parameter                  :: pi=3.1415926535897932
 character*1 :: tilef

! Remove any halo region.

 i_target_out = i_target-(2*halo)
 j_target_out = j_target-(2*halo)
 
 i_start = halo + 1
 j_start = halo + 1
 i_end   = i_target - halo
 j_end   = j_target - halo

 ip1_target_out = i_target_out + 1
 jp1_target_out = j_target_out + 1

 ip1_end = i_end + 1
 jp1_end = j_end + 1

 nlev=lev_target
 print*," SG SG SG NUMLEVELS=",nlev



 allocate(x_data(i_target_out))
 do i = 1, i_target_out
   x_data(i) = float(i)
 enddo

 allocate(x_datap1(ip1_target_out))
 do i = 1, ip1_target_out
   x_datap1(i) = float(i)
 enddo

 allocate(y_data(j_target_out))
 do i = 1, j_target_out
   y_data(i) = float(i)
 enddo

 allocate(y_datap1(jp1_target_out))
 do i = 1, jp1_target_out
   y_datap1(i) = float(i)
 enddo

 allocate(lev_data(lev_target))
 do i = 1, lev_target
   lev_data(i) = float(i)
 enddo
 print*,"LEVDATA ", lev_data


 if (localpet == 0) then
   allocate(data_one_tile(i_target,j_target))
   allocate(data_one_tile_ijp1(i_target,jp1_target))
   allocate(data_one_tile_ip1j(ip1_target,j_target))
   allocate(data_one_tile_ip1jp1(ip1_target,jp1_target))
   allocate(dum2d_ijp1(i_target_out,jp1_target_out))
   allocate(dum2d_ip1j(ip1_target_out,j_target_out))
   allocate(dum2d_ip1jp1(ip1_target_out,jp1_target_out))
   allocate(dum2d(i_target_out,j_target_out))
 else
   allocate(data_one_tile(0,0))
   allocate(data_one_tile_ijp1(0,0))
   allocate(data_one_tile_ip1j(0,0))
   allocate(data_one_tile_ip1jp1(0,0))
   allocate(dum2d_ijp1(0,0))
   allocate(dum2d_ip1j(0,0))
   allocate(dum2d_ip1jp1(0,0))
   allocate(dum2d(0,0))
 endif



 TILE_LOOP : do tile = 1, num_tiles_target_grid

   LOCAL_PET : if (localpet == 0) then

!     WRITE(OUTFILE, '(A, I1, A)'), 'out.atm.tile', tile, '.nc'
     WRITE(OUTFILE, '(A, I1, A)'), 'fv_core.res.tile', tile, '.nc'

!--- open the file
     error = nf90_create(outfile, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                         ncid, initialsize=inital, chunksize=fsize)
     call netcdf_err(error, 'CREATING FILE='//trim(outfile) )

!--- define dimension
     error = nf90_def_dim(ncid, 'xaxis_1', i_target_out, dim_lon)
     call netcdf_err(error, 'DEFINING LON DIMENSION' )
     error = nf90_def_dim(ncid, 'xaxis_2', ip1_target_out, dim_lonp)
     call netcdf_err(error, 'DEFINING LONP DIMENSION' )
     error = nf90_def_dim(ncid, 'yaxis_1', jp1_target_out, dim_latp)
     call netcdf_err(error, 'DEFINING LATP DIMENSION' )
     error = nf90_def_dim(ncid, 'yaxis_2', j_target_out, dim_lat)
     call netcdf_err(error, 'DEFINING LAT DIMENSION' )
     error = nf90_def_dim(ncid, 'zaxis_1', lev_target, dim_lev)
     call netcdf_err(error, 'DEFINING LEV DIMENSION' )
     error = nf90_def_dim(ncid, 'Time', 1, dim_time)
     call netcdf_err(error, 'DEFINING TIME DIMENSION' )

 !--- define field
     error = nf90_def_var(ncid, 'xaxis_1', NF90_DOUBLE, (/dim_lon/), id_x)
     call netcdf_err(error, 'DEFINING x1 LON FIELD' )
     error = nf90_put_att(ncid, id_x, "long_name", "xaxis_1")
     error = nf90_put_att(ncid, id_x, "units", "none")
     error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
     call netcdf_err(error, 'WRITING x1 LON FIELD' )

     error = nf90_def_var(ncid, 'xaxis_2', NF90_DOUBLE, (/dim_lonp/), id_xp)
     call netcdf_err(error, 'DEFINING x2 LON FIELD' )
     error = nf90_put_att(ncid, id_xp, "long_name", "xaxis_2")
     error = nf90_put_att(ncid, id_xp, "units", "none")
     error = nf90_put_att(ncid, id_xp, "cartesian_axis", "X")
     call netcdf_err(error, 'WRITING x2 LON FIELD' )

     error = nf90_def_var(ncid, 'yaxis_1', NF90_DOUBLE, (/dim_latp/), id_yp)
     call netcdf_err(error, 'DEFINING y1 LAT FIELD' )
     error = nf90_put_att(ncid, id_yp, "long_name", "yaxis_1")
     error = nf90_put_att(ncid, id_yp, "units", "none")
     error = nf90_put_att(ncid, id_yp, "cartesian_axis", "Y")
     call netcdf_err(error, 'WRITING y1 LAT FIELD' )

     error = nf90_def_var(ncid, 'yaxis_2', NF90_DOUBLE, (/dim_lat/), id_y)
     call netcdf_err(error, 'DEFINING y2 LAT FIELD' )
     error = nf90_put_att(ncid, id_y, "long_name", "yaxis_2")
     error = nf90_put_att(ncid, id_y, "units", "none")
     error = nf90_put_att(ncid, id_y, "cartesian_axis", "Y")
     call netcdf_err(error, 'WRITING y2 LAT FIELD' )

     error = nf90_def_var(ncid, 'zaxis_1', NF90_DOUBLE, (/dim_lev/), id_z)
     call netcdf_err(error, 'DEFINING z1 lev FIELD' )
     error = nf90_put_att(ncid, id_z, "long_name", "zaxis_1")
     error = nf90_put_att(ncid, id_z, "units", "none")
     error = nf90_put_att(ncid, id_z, "cartesian_axis", "Z")
     call netcdf_err(error, 'WRITING z1 lev FIELD' )

     error = nf90_def_var(ncid, 'Time', NF90_DOUBLE, dim_time, id_time)
     call netcdf_err(error, 'DEFINING TIME FIELD' )
     error = nf90_put_att(ncid, id_time, "long_name", "Time")
     error = nf90_put_att(ncid, id_time, "units", "time level")
     error = nf90_put_att(ncid, id_time, "cartesian_axis", "T")
     call netcdf_err(error, 'WRITING TIME FIELD' )



     error = nf90_def_var(ncid, 'u', NF90_FLOAT,(/dim_lon,dim_latp,dim_lev,dim_time/),id_u)
     error = nf90_put_att(ncid, id_u, "long_name", "u wind")
     error = nf90_put_att(ncid, id_u, "units", "m/s")
     call netcdf_err(error, 'WRITING U' )

     error = nf90_def_var(ncid, 'v', NF90_FLOAT,(/dim_lonp,dim_lat,dim_lev,dim_time/),id_v)
     error = nf90_put_att(ncid, id_v, "long_name", "v wind")
     error = nf90_put_att(ncid, id_v, "units", "m/s")
     call netcdf_err(error, 'WRITING V' )





!     error = nf90_def_var(ncid, 'u', NF90_DOUBLE, (/dim_lon,dim_latp,dim_lev,dim_time/),id_u_s)
!     call netcdf_err(error, 'WRITING U_S' )
!     error = nf90_put_att(ncid, id_u_s, "long_name", "u")
!     error = nf90_put_att(ncid, id_u_s, "units", "none")
!     call netcdf_err(error, 'WRITING U_S' )
!
!     error = nf90_def_var(ncid, 'v', NF90_DOUBLE, (/dim_lonp,dim_lat,dim_lev,dim_time/),id_v_w)
!     call netcdf_err(error, 'WRITING V_W' )
!     error = nf90_put_att(ncid, id_v_w, "long_name", "v")
!     error = nf90_put_att(ncid, id_v_w, "units", "none")
!     call netcdf_err(error, 'WRITING V_W' )

     error = nf90_def_var(ncid, 'W', NF90_DOUBLE, (/dim_lon,dim_lat,dim_lev,dim_time/), id_w)
     call netcdf_err(error, 'WRITING W' )
     error = nf90_put_att(ncid, id_w, "long_name", "W")
     error = nf90_put_att(ncid, id_w, "units", "none")
     call netcdf_err(error, 'WRITING W' )

     error = nf90_def_var(ncid, 'T', NF90_DOUBLE, (/dim_lon,dim_lat,dim_lev,dim_time/), id_t)
     call netcdf_err(error, 'WRITING T' )
     error = nf90_put_att(ncid, id_t, "long_name", "T")
     error = nf90_put_att(ncid, id_t, "units", "none")
     call netcdf_err(error, 'WRITING T' )

     error = nf90_def_var(ncid, 'delp', NF90_DOUBLE, (/dim_lon,dim_lat,dim_lev,dim_time/), id_delp)
     call netcdf_err(error, 'WRITING DELP' )
     error = nf90_put_att(ncid, id_delp, "long_name", "delp")
     error = nf90_put_att(ncid, id_delp, "units", "none")
     call netcdf_err(error, 'WRITING T' )

     error = nf90_def_var(ncid, 'phis', NF90_DOUBLE, (/dim_lon,dim_lat,dim_time/), id_phis)
     call netcdf_err(error, 'WRITING phis' )
     error = nf90_put_att(ncid, id_phis, "long_name", "phis")
     error = nf90_put_att(ncid, id_phis, "units", "none")
     call netcdf_err(error, 'WRITING phis' )


     error = nf90_enddef(ncid, header_buffer_val,4,0,4)
     call netcdf_err(error, 'DEFINING HEADER' )

   endif LOCAL_PET ! is localpet 0?

!------------------------------------ 

   if (localpet == 0) then
     error = nf90_put_var( ncid, id_z, lev_data)
     call netcdf_err(error, 'WRITING ZAXIS RECORD' )
     error = nf90_put_var( ncid, id_x, x_data)
     call netcdf_err(error, 'WRITING XAXIS RECORD' )
     error = nf90_put_var( ncid, id_xp, x_datap1)
     call netcdf_err(error, 'WRITING XAXISp RECORD' )
     error = nf90_put_var( ncid, id_y, y_data)
     call netcdf_err(error, 'WRITING YAXIS RECORD' )
     error = nf90_put_var( ncid, id_yp, y_datap1)
     call netcdf_err(error, 'WRITING YAXISp RECORD' )
     times = 1.0
     error = nf90_put_var( ncid, id_time, times)
     call netcdf_err(error, 'WRITING TIME RECORD' )
   endif

!  phis

   print*,"- CALL FieldGather FOR TERRAIN TARGET GRID to calculate PHIS: ", tile
   call ESMF_FieldGather(terrain_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(i_start:i_end, j_start:j_end)
     error = nf90_put_var( ncid, id_phis, gravity*dum2d)
     call netcdf_err(error, 'WRITING PHIS RECORD' )
   endif

!  vertical velocity

   if (localpet == 0) then
     allocate(dum3d(i_target_out,j_target_out,lev_target))
     allocate(data_one_tile_3d(i_target,j_target,lev_target))
   else
     allocate(dum3d(0,0,0))
     allocate(data_one_tile_3d(0,0,0))
   endif

   print*,"- CALL FieldGather FOR TARGET GRID VERTICAL VELOCITY FOR TILE: ", tile
   call ESMF_FieldGather(dzdt_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
     error = nf90_put_var( ncid, id_w, dum3d)
     call netcdf_err(error, 'WRITING VERTICAL VELOCITY RECORD' )
   endif

!  delp

   print*,"- CALL FieldGather FOR TARGET GRID DELP FOR TILE: ", tile
   call ESMF_FieldGather(delp_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
     error = nf90_put_var( ncid, id_delp, dum3d)
     call netcdf_err(error, 'WRITING DELP RECORD' )
   endif


!  temperature

   print*,"- CALL FieldGather FOR TARGET GRID TEMPERATURE FOR TILE: ", tile
   call ESMF_FieldGather(temp_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
     error = nf90_put_var( ncid, id_t, dum3d)
     call netcdf_err(error, 'WRITING TEMPERTAURE RECORD' )
   endif


   deallocate(dum3d, data_one_tile_3d)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  latlon arrays for rotating uv
   print*,"localpet=",localpet
   if (localpet == 0) then
     allocate(lat_corner(ip1_target,jp1_target))
     allocate(lon_corner(ip1_target,jp1_target))
   else
     allocate(lat_corner(0,0))
     allocate(lon_corner(0,0))
   endif

!  longitude_corners

   print*,"- CALL FieldGather FOR TARGET GRID LONGITUDE_c FOR TILE: ", tile
   call ESMF_FieldGather(longitude_c_target_grid, data_one_tile_ip1jp1, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d_ip1jp1(:,:) = data_one_tile_ip1jp1(i_start:ip1_end, j_start:jp1_end)
     lon_corner=dum2d_ip1jp1*pi/180
   endif
   print*,"SG SG SG SG shape(lonCORNER):",shape(lon_corner)
   print*,"SG SG SG lonCORNER=",lon_corner

!  latitude_corners

   print*,"- CALL FieldGather FOR TARGET GRID LATITUDE_c FOR TILE: ", tile
   call ESMF_FieldGather(latitude_c_target_grid, data_one_tile_ip1jp1, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d_ip1jp1(:,:) = data_one_tile_ip1jp1(i_start:ip1_end, j_start:jp1_end)
     lat_corner=dum2d_ip1jp1*pi/180
   endif
   print*,"SG SG SG SG shape(latCORNER):",shape(lat_corner)
   print*,"SG SG SG latCORNER=",lat_corner

!  uwinds s

   if (localpet == 0) then
     allocate(dum3d(i_target_out,jp1_target_out,lev_target))
     allocate(u_s(i_target_out,jp1_target_out,lev_target))
     allocate(v_s(i_target_out,jp1_target_out,lev_target))
     allocate(data_one_tile_3d(i_target,jp1_target,lev_target))
   else
     allocate(dum3d(0,0,0))
     allocate(u_s(0,0,0))
     allocate(v_s(0,0,0))
     allocate(data_one_tile_3d(0,0,0))
   endif

   print*,"- CALL FieldGather FOR TARGET GRID U_S FOR TILE: ", tile
   call ESMF_FieldGather(u_s_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, jp1_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:jp1_end,:)
     u_s=dum3d(:,:,:)
   endif

!  vwinds s

   print*,"- CALL FieldGather FOR TARGET GRID V_S FOR TILE: ", tile
   call ESMF_FieldGather(v_s_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, jp1_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:jp1_end,:)
     v_s=dum3d(:,:,:)
   endif

   deallocate(dum3d, data_one_tile_3d)


!  uwinds w

   if (localpet == 0) then
     allocate(dum3d(ip1_target_out,j_target_out,lev_target))
     allocate(u_w(ip1_target_out,j_target_out,lev_target))
     allocate(v_w(ip1_target_out,j_target_out,lev_target))
     allocate(data_one_tile_3d(ip1_target,j_target,lev_target))
   else
     allocate(dum3d(0,0,0))
     allocate(u_w(0,0,0))
     allocate(v_w(0,0,0))
     allocate(data_one_tile_3d(0,0,0))
   endif

   print*,"- CALL FieldGather FOR TARGET GRID U_W FOR TILE: ", tile
   call ESMF_FieldGather(u_w_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, ip1_target, j_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:ip1_end,j_start:j_end,:)
     u_w=dum3d(:,:,:)
   endif

!  vwinds w

   print*,"- CALL FieldGather FOR TARGET GRID V_W FOR TILE: ", tile
   call ESMF_FieldGather(v_w_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, ip1_target, j_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:ip1_end,j_start:j_end,:)
     v_w=dum3d(:,:,:)
   endif

   print*,"- SG SG SG SG SG got here "
   deallocate(dum3d, data_one_tile_3d)
   print*,"- SG SG SG SG SG got here 2"

!-------------------------------------------------------------------------------
! v point (left side)
! vector stuff

!   nlev=shape(u_w,3)
   print*," SG SG SG NUMLEVELS=",nlev
   if (localpet == 0) then
     allocate(p1(2))
     allocate(p2(2))
     allocate(p3(2))
     print*,"- SG SG SG SG SG p1p2p3: "
     allocate(e1(3))
     allocate(ex(3))
     allocate(ey(3))
     print*,"- SG SG SG SG SG e1exey: "
     allocate(udnew(i_target,j_target+1,nlev))
     allocate(vdnew(i_target+1,j_target ,nlev ))
     print*,"- SG SG SG before vdnew"

     do j=1,j_target
        do i=1,i_target+1
           p1(1) = lon_corner(i,j)
           p1(2) = lat_corner(i,j)
           p2(1) = lon_corner(i,j+1)
           p2(2) = lat_corner(i,j+1)
           call mid_pt_sphere(p1, p2, p3)
           call get_unit_vect2(p1, p2, e1)
           call get_latlon_vector(p3, ex, ey)
           print*,"- SG SG SG SG ijk: ", shape(u_w)
           print*,"SG SG SG SG dimlev=",nlev
           do k=1,nlev
              vdnew(i,j,k) = u_w(i,j,k)*inner_prod(e1,ex) + &
                             v_w(i,j,k)*inner_prod(e1,ey)
           enddo
        enddo
     enddo
! upoint south
     print*,"- SG SG SG between vdnew udnew"
     do j=1,j_target+1
        do i=1,i_target
           p1(1) = lon_corner(i,j)
           p1(2) = lat_corner(i,j)
           p2(1) = lon_corner(i+1,j)
           p2(2) = lat_corner(i+1,j)
           call mid_pt_sphere(p1, p2, p3)
           call get_unit_vect2(p1, p2, e1)
           call get_latlon_vector(p3, ex, ey)
           do k=1,nlev
              udnew(i,j,k) = u_s(i,j,k)*inner_prod(e1,ex) + &
                             v_s(i,j,k)*inner_prod(e1,ey)
           enddo
        enddo
     enddo
     print*,"- SG SG SG SG SG UVdnew/UV_wsavg=",sqrt(udnew**2+vdnew**2)/((sqrt(u_s**2 + v_s**2)+sqrt(u_w**2 + v_w**2))/2)
   else
     allocate(p1(0))
     allocate(p2(0))
     allocate(p3(0))
     allocate(e1(0))
     allocate(ex(0))
     allocate(ey(0))
     allocate(udnew(0,0,0))
     allocate(vdnew(0,0,0))

     print*,"- SG SG SG after udnew"
   endif
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! WRITE U

   if (localpet == 0) then
     allocate(dum3d(i_target_out,jp1_target_out,lev_target))
     allocate(data_one_tile_3d(i_target,jp1_target,lev_target))
   else
     allocate(dum3d(0,0,0))
     allocate(data_one_tile_3d(0,0,0))
   endif

   print*,"- making NEW U: ", tile

   if (localpet == 0) then
     error = nf90_put_var( ncid, id_u, udnew)
     call netcdf_err(error, 'WRITING U_new RECORD' )
   endif

   deallocate(dum3d, data_one_tile_3d)

! WRITE V 

   if (localpet == 0) then
     allocate(dum3d(ip1_target_out,j_target_out,lev_target))
     allocate(data_one_tile_3d(ip1_target,j_target,lev_target))
   else
     allocate(dum3d(0,0,0))
     allocate(data_one_tile_3d(0,0,0))
   endif

   print*,"- making NEW V: ", tile

   if (localpet == 0) then
     error = nf90_put_var( ncid, id_v, vdnew)
     call netcdf_err(error, 'WRITING V_new RECORD' )
   endif

   deallocate(dum3d, data_one_tile_3d)

   print*,"- GOT after NEW V ", tile

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   deallocate(u_s, u_w, v_s, v_w)
   deallocate(p1,p2,p3,e1,ex,ey)
   deallocate(udnew,vdnew)
   deallocate(lat_corner,lon_corner)

!-------------------------------------------------------------------------------
   print*,"- GOT to end of core ", tile

! close file
!-------------------------------------------------------------------------------

  if (localpet == 0) error = nf90_close(ncid)
  print*,"- closed core ", tile

!-------------------------------------------------------------------------------

   LOCAL_PET2 : if (localpet == 0) then

     WRITE(OUTFILE, '(A, I1, A)'), 'fv_tracer.res.tile', tile, '.nc'


     print*,"- GOT TO TRACER ", tile


!--- open the file
     error = nf90_create(outfile, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                         ncid, initialsize=inital, chunksize=fsize)
     call netcdf_err(error, 'CREATING FILE='//trim(outfile) )


!--- define dimension
     error = nf90_def_dim(ncid, 'xaxis_1', i_target_out, dim_lon)
     call netcdf_err(error, 'DEFINING x1 DIMENSION' )
     error = nf90_def_dim(ncid, 'yaxis_1', j_target_out, dim_lat)
     call netcdf_err(error, 'DEFINING y1 DIMENSION' )
     error = nf90_def_dim(ncid, 'zaxis_1', lev_target, dim_lev)
     call netcdf_err(error, 'DEFINING LEV DIMENSION' )
     error = nf90_def_dim(ncid, 'Time', 1, dim_time)
     call netcdf_err(error, 'DEFINING TIME DIMENSION' )

 !--- define field
     error = nf90_def_var(ncid, 'xaxis_1', NF90_DOUBLE, (/dim_lon/), id_x)
     call netcdf_err(error, 'DEFINING x1 LON FIELD' )
     error = nf90_put_att(ncid, id_x, "long_name", "xaxis_1")
     error = nf90_put_att(ncid, id_x, "units", "none")
     error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
     call netcdf_err(error, 'WRITING x1 LON FIELD' )

     error = nf90_def_var(ncid, 'yaxis_1', NF90_DOUBLE, (/dim_lat/), id_y)
     call netcdf_err(error, 'DEFINING y LAT FIELD' )
     error = nf90_put_att(ncid, id_y, "long_name", "yaxis_1")
     error = nf90_put_att(ncid, id_y, "units", "none")
     error = nf90_put_att(ncid, id_y, "cartesian_axis", "Y")
     call netcdf_err(error, 'WRITING y1 LAT FIELD' )

     error = nf90_def_var(ncid, 'zaxis_1', NF90_DOUBLE, (/dim_lev/), id_z)
     call netcdf_err(error, 'DEFINING z1 lev FIELD' )
     error = nf90_put_att(ncid, id_z, "long_name", "zaxis_1")
     error = nf90_put_att(ncid, id_z, "units", "none")
     error = nf90_put_att(ncid, id_z, "cartesian_axis", "Z")
     call netcdf_err(error, 'WRITING z1 lev FIELD' )

     error = nf90_def_var(ncid, 'Time', NF90_DOUBLE, dim_time, id_time)
     call netcdf_err(error, 'DEFINING TIME FIELD' )
     error = nf90_put_att(ncid, id_time, "long_name", "Time")
     error = nf90_put_att(ncid, id_time, "units", "time level")
     error = nf90_put_att(ncid, id_time, "cartesian_axis", "T")
     call netcdf_err(error, 'WRITING TIME FIELD' )


     error = nf90_def_var(ncid, 'sphum', NF90_FLOAT,(/dim_lon,dim_lat,dim_lev,dim_time/), id_sphum)
     call netcdf_err(error, 'WRITING SPHUM' )
     error = nf90_put_att(ncid, id_sphum, "long_name", "sphum")
     error = nf90_put_att(ncid, id_sphum, "units", "none")
     call netcdf_err(error, 'WRITING SPHUM' )

     error = nf90_def_var(ncid, 'liq_wat',NF90_FLOAT,(/dim_lon,dim_lat,dim_lev,dim_time/), id_liq_wat)
     call netcdf_err(error, 'WRITING liq_wat' )
     error = nf90_put_att(ncid, id_liq_wat, "long_name", "liq_wat")
     error = nf90_put_att(ncid, id_liq_wat, "units", "none")
     call netcdf_err(error, 'WRITING liq_wat' )


     if (gfdl_mp) then
       error = nf90_def_var(ncid, 'rainwat',NF90_FLOAT,(/dim_lon,dim_lat,dim_lev,dim_time/), id_rwmr)
       call netcdf_err(error, 'WRITING rainwat' )
       error = nf90_put_att(ncid, id_rwmr, "long_name", "rainwat")
       error = nf90_put_att(ncid, id_rwmr, "units", "none")
       call netcdf_err(error, 'WRITING rainwat' )

       error = nf90_def_var(ncid, 'ice_wat',NF90_FLOAT,(/dim_lon,dim_lat,dim_lev,dim_time/), id_icmr)
       call netcdf_err(error, 'WRITING ice_wat' )
       error = nf90_put_att(ncid, id_icmr, "long_name", "ice_wat")
       error = nf90_put_att(ncid, id_icmr, "units", "none")
       call netcdf_err(error, 'WRITING ice_wat' )

       error = nf90_def_var(ncid, 'snowwat',NF90_FLOAT,(/dim_lon,dim_lat,dim_lev,dim_time/), id_snmr)
       call netcdf_err(error, 'WRITING snowwat' )
       error = nf90_put_att(ncid, id_snmr, "long_name", "snowwat")
       error = nf90_put_att(ncid, id_snmr, "units", "none")
       call netcdf_err(error, 'WRITING snowwat' )

       error = nf90_def_var(ncid, 'graupel',NF90_FLOAT,(/dim_lon,dim_lat,dim_lev,dim_time/), id_grle)
       call netcdf_err(error, 'WRITING graupel' )
       error = nf90_put_att(ncid, id_grle, "long_name", "graupel")
       error = nf90_put_att(ncid, id_grle, "units", "none")
       call netcdf_err(error, 'WRITING graupel' )
     endif

     error = nf90_def_var(ncid, 'o3mr', NF90_FLOAT,(/dim_lon,dim_lat,dim_lev,dim_time/),id_o3mr)
     call netcdf_err(error, 'WRITING o3mr' )
     error = nf90_put_att(ncid, id_o3mr, "long_name", "o3mr")
     error = nf90_put_att(ncid, id_o3mr, "units", "none")
     call netcdf_err(error, 'WRITING o3mr' )


     error = nf90_enddef(ncid, header_buffer_val,4,0,4)
     call netcdf_err(error, 'DEFINING HEADER' )

   endif LOCAL_PET2 ! is localpet 0?
!---
   if (localpet == 0) then
     error = nf90_put_var( ncid, id_z, lev_data)
     call netcdf_err(error, 'WRITING ZAXIS RECORD' )
     error = nf90_put_var( ncid, id_x, x_data)
     call netcdf_err(error, 'WRITING XAXIS RECORD' )
     error = nf90_put_var( ncid, id_y, y_data)
     call netcdf_err(error, 'WRITING YAXIS RECORD' )
     times = 1.0
     error = nf90_put_var( ncid, id_time, times)
     call netcdf_err(error, 'WRITING TIME RECORD' )
   endif


   if (localpet == 0) then
     allocate(dum3d(i_target_out,j_target_out,lev_target))
     allocate(data_one_tile_3d(i_target,j_target,lev_target))
   else
     allocate(dum3d(0,0,0))
     allocate(data_one_tile_3d(0,0,0))
   endif

!  specific humidity

   print*,"- CALL FieldGather FOR TARGET GRID SPECIFIC HUMIDITY FOR TILE: ",tile
   call ESMF_FieldGather(spec_hum_target_grid, data_one_tile_3d,rootPet=0,tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
     error = nf90_put_var( ncid, id_sphum, dum3d)
     call netcdf_err(error, 'WRITING SPECIFIC HUMIDITY RECORD' )
   endif

!  ozone

   print*,"- CALL FieldGather FOR TARGET GRID OZONE MIXING RATIO FOR TILE:",tile
   call ESMF_FieldGather(o3mr_target_grid, data_one_tile_3d,rootPet=0,tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
     error = nf90_put_var( ncid, id_o3mr, dum3d)
     call netcdf_err(error, 'WRITING OZONE MIXING RATIO RECORD' )
   endif

!  liquid water

   print*,"- CALL FieldGather FOR TARGET GRID LIQUID WATER MIXING RATIO FOR TILE: ", tile
   call ESMF_FieldGather(liq_wat_target_grid, data_one_tile_3d,rootPet=0,tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
      call mpi_abort

   if (localpet == 0) then
     call flip(data_one_tile_3d, i_target, j_target, lev_target)
     dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
     error = nf90_put_var( ncid, id_liq_wat, dum3d)
     call netcdf_err(error, 'WRITING LIQUID WATER MIXING RATIO RECORD' )
   endif
   if (gfdl_mp) then

!    rain water

     print*,"- CALL FieldGather FOR TARGET GRID RAIN WATER MIXING RATIO FOR TILE: ", tile
     call ESMF_FieldGather(rwmr_target_grid, data_one_tile_3d,rootPet=0,tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
       call mpi_abort

     if (localpet == 0) then
       call flip(data_one_tile_3d, i_target, j_target, lev_target)
       dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
       error = nf90_put_var( ncid, id_rwmr, dum3d)
       call netcdf_err(error, 'WRITING RAIN WATER MIXING RATIO RECORD' )
     endif

!    ice water

     print*,"- CALL FieldGather FOR TARGET GRID ICE WATER MIXING RATIO FOR TILE:", tile
     call ESMF_FieldGather(icmr_target_grid, data_one_tile_3d,rootPet=0,tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
       call mpi_abort

     if (localpet == 0) then
       call flip(data_one_tile_3d, i_target, j_target, lev_target)
       dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
       error = nf90_put_var( ncid, id_icmr, dum3d)
       call netcdf_err(error, 'WRITING ICE WATER MIXING RATIO RECORD' )
     endif

!    snow water

     print*,"- CALL FieldGather FOR TARGET GRID SNOW WATER MIXING RATIO FOR TILE: ", tile
     call ESMF_FieldGather(snmr_target_grid, data_one_tile_3d,rootPet=0,tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
       call mpi_abort

     if (localpet == 0) then
       call flip(data_one_tile_3d, i_target, j_target, lev_target)
       dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
       error = nf90_put_var( ncid, id_snmr, dum3d)
       call netcdf_err(error, 'WRITING SNOW WATER MIXING RATIO RECORD' )
     endif

!    graupel

     print*,"- CALL FieldGather FOR TARGET GRID GRAUPEL FOR TILE: ", tile
     call ESMF_FieldGather(grle_target_grid, data_one_tile_3d,rootPet=0,tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
       call mpi_abort

     if (localpet == 0) then
       call flip(data_one_tile_3d, i_target, j_target, lev_target)
       dum3d(:,:,:) = data_one_tile_3d(i_start:i_end,j_start:j_end,:)
       error = nf90_put_var( ncid, id_grle, dum3d)
       call netcdf_err(error, 'WRITING GRAUPEL RECORD' )
     endif

   endif

   deallocate(dum3d, data_one_tile_3d)

!-------------------------------------------------------------------------------
! close file
!-------------------------------------------------------------------------------

   if (localpet == 0) error = nf90_close(ncid)


 enddo TILE_LOOP

 deallocate(data_one_tile)
 deallocate(dum2d)

 end subroutine write_fv3_atm_data_netcdf


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


subroutine mid_pt_sphere(p1, p2, pm)
implicit none
real , intent(IN)  :: p1(2), p2(2)
real , intent(OUT) :: pm(2)
real e1(3), e2(3), e3(3)

call latlon2xyz(p1, e1)
call latlon2xyz(p2, e2)
call mid_pt3_cart(e1, e2, e3)
call cart_to_latlon(1, e3, pm(1), pm(2))

end subroutine mid_pt_sphere

subroutine get_unit_vect2( e1, e2, uc )
implicit none
real, intent(in) :: e1(2), e2(2)
real, intent(out):: uc(3) !< unit vector e1--->e2
! Local:
real, dimension(3):: pc, p1, p2, p3

call latlon2xyz(e1, p1)
call latlon2xyz(e2, p2)

call mid_pt3_cart(p1, p2,  pc)
call vect_cross(p3, p2, p1)
call vect_cross(uc, pc, p3)
call normalize_vect( uc )

end subroutine get_unit_vect2

subroutine mid_pt3_cart(p1, p2, e)
implicit none
real, intent(IN)  :: p1(3), p2(3)
real, intent(OUT) :: e(3)
!
real :: q1(3), q2(3)
real :: dd, e1, e2, e3
integer k

do k=1,3
   q1(k) = p1(k)
   q2(k) = p2(k)
enddo

e1 = q1(1) + q2(1)
e2 = q1(2) + q2(2)
e3 = q1(3) + q2(3)

dd = sqrt( e1**2 + e2**2 + e3**2 )
e1 = e1 / dd
e2 = e2 / dd
e3 = e3 / dd

e(1) = e1
e(2) = e2
e(3) = e3

end subroutine mid_pt3_cart


subroutine get_latlon_vector(pp, elon, elat)
implicit none
real, intent(IN)  :: pp(2)
real, intent(OUT) :: elon(3), elat(3)

elon(1) = -SIN(pp(1))
elon(2) =  COS(pp(1))
elon(3) =  0.0
elat(1) = -SIN(pp(2))*COS(pp(1))
elat(2) = -SIN(pp(2))*SIN(pp(1))
elat(3) =  COS(pp(2))

end subroutine get_latlon_vector
subroutine latlon2xyz(p, e)
implicit none
real, intent(in) :: p(2)
real, intent(out):: e(3)


e(1) = cos(p(2)) * cos(p(1))
e(2) = cos(p(2)) * sin(p(1))
e(3) = sin(p(2))

end subroutine latlon2xyz



!>@brief The subroutine 'vect_cross performs cross products
!! of 3D vectors: e = P1 X P2
subroutine vect_cross(e, p1, p2)
implicit none
real, intent(in) :: p1(3), p2(3)
real, intent(out):: e(3)

e(1) = p1(2)*p2(3) - p1(3)*p2(2)
e(2) = p1(3)*p2(1) - p1(1)*p2(3)
e(3) = p1(1)*p2(2) - p1(2)*p2(1)

end subroutine vect_cross


!>@brief The subroutine 'normalize_vect' makes 'e' a unit vector.
subroutine normalize_vect(e)
implicit none
real, intent(inout):: e(3)
real:: pdot
integer k

pdot = e(1)**2 + e(2)**2 + e(3)**2
pdot = sqrt( pdot )

do k=1,3
   e(k) = e(k) / pdot
enddo

end subroutine normalize_vect


subroutine cart_to_latlon(np, q, xs, ys)
implicit none
! vector version of cart_to_latlon1
integer, intent(in):: np
real, intent(inout):: q(3,np)
real, intent(inout):: xs(np), ys(np)
! local
real, parameter:: esl=1.d-10
real, parameter:: pi=3.1415926535897932
real :: p(3)
real :: dist, lat, lon
integer i,k

do i=1,np
   do k=1,3
      p(k) = q(k,i)
   enddo
   dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
   do k=1,3
      p(k) = p(k) / dist
   enddo
   if ( (abs(p(1))+abs(p(2)))  < esl ) then
        lon = 0.0
   else
        lon = atan2( p(2), p(1) )   ! range [-pi,pi]
   endif

   if ( lon < 0.) lon = 2.0*pi + lon

   lat = asin(p(3))

   xs(i) = lon
   ys(i) = lat

! q Normalized:
   do k=1,3
      q(k,i) = p(k)
   enddo
enddo

end  subroutine cart_to_latlon

real function inner_prod(v1, v2)
implicit none
real,intent(in):: v1(3), v2(3)
real :: prod16
integer k

prod16 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
inner_prod = prod16

end function inner_prod
 


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

 subroutine write_fv3_sfc_data_netcdf(localpet)

 use esmf
 use netcdf

 use model_grid, only            : num_tiles_target_grid, &
                                   landmask_target_grid, &
                                   i_target, j_target, lsoil_target

 use program_setup, only         : convert_nst, halo

 use surface, only               : canopy_mc_target_grid,  &
                                   f10m_target_grid, &
                                   ffmm_target_grid, &
                                   q2m_target_grid,   &
                                   seaice_depth_target_grid, &
                                   seaice_fract_target_grid, &
                                   seaice_skin_temp_target_grid, &
                                   skin_temp_target_grid, &
                                   soil_temp_target_grid, &
                                   soilm_liq_target_grid, &
                                   soilm_tot_target_grid, &
                                   srflag_target_grid, &
                                   snow_liq_equiv_target_grid, &
                                   snow_depth_target_grid, &
                                   t2m_target_grid,   &
                                   tprcp_target_grid, &
                                   ustar_target_grid, &
                                   z0_target_grid, &
                                   c_d_target_grid, &
                                   c_0_target_grid, &
                                   d_conv_target_grid, &
                                   dt_cool_target_grid, &
                                   ifd_target_grid, &
                                   qrain_target_grid, &
                                   tref_target_grid, &
                                   w_d_target_grid, &
                                   w_0_target_grid, &
                                   xs_target_grid, &
                                   xt_target_grid, &
                                   xu_target_grid, &
                                   xv_target_grid, &
                                   xz_target_grid, &
                                   xtts_target_grid, &
                                   xzts_target_grid, &
                                   z_c_target_grid, &
                                   zm_target_grid

 use static_data, only           : alvsf_target_grid,   &
                                   alnsf_target_grid,   &
                                   alvwf_target_grid,   &
                                   alnwf_target_grid,   &
                                   facsf_target_grid, &
                                   facwf_target_grid, &
                                   max_veg_greenness_target_grid, &
                                   min_veg_greenness_target_grid, &
                                   mxsno_albedo_target_grid, &
                                   slope_type_target_grid, &
                                   soil_type_target_grid,  &
                                   substrate_temp_target_grid,  &
                                   veg_greenness_target_grid, &
                                   veg_type_target_grid

 implicit none

 integer, intent(in)            :: localpet
 character(len=128)             :: outfile

 integer                        :: fsize=65536, inital = 0
 integer                        :: header_buffer_val = 16384
 integer                        :: dim_x, dim_y, dim_lsoil, dim_time
 integer                        :: error, i, ncid, tile
 integer                        :: id_x, id_y, id_lsoil
 integer                        :: id_slmsk, id_time
 integer                        :: id_tsea, id_sheleg, id_tg3
 integer                        :: id_zorl, id_alvsf, id_alvwf
 integer                        :: id_alnsf, id_alnwf, id_vfrac
 integer                        :: id_canopy, id_f10m, id_t2m
 integer                        :: id_q2m, id_vtype, id_stype
 integer                        :: id_facsf, id_facwf, id_uustar
 integer                        :: id_ffmm, id_ffhh, id_hice
 integer                        :: id_fice, id_tisfc, id_tprcp
 integer                        :: id_srflag, id_snwdph, id_shdmin
 integer                        :: id_shdmax, id_slope, id_snoalb
 integer                        :: id_stc, id_smc, id_slc
 integer                        :: id_tref, id_z_c, id_c_0
 integer                        :: id_c_d, id_w_0, id_w_d
 integer                        :: id_xt, id_xs, id_xu, id_xv
 integer                        :: id_xz, id_zm, id_xtts, id_xzts
 integer                        :: id_d_conv, id_ifd, id_dt_cool
 integer                        :: id_qrain
 integer                        :: i_target_out, j_target_out
 integer                        :: istart, iend, jstart, jend

 integer(esmf_kind_i8), allocatable :: idata_one_tile(:,:)

 real(kind=4), allocatable       :: lsoil_data(:), x_data(:), y_data(:)
 real(kind=8), allocatable       :: dum2d(:,:), dum3d(:,:,:)
 real(kind=4)                    :: times
 real(esmf_kind_r8), allocatable :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable :: data_one_tile_3d(:,:,:)

! Remove any halo region.

 i_target_out = i_target-(2*halo)
 j_target_out = j_target-(2*halo)
 
 istart = halo + 1
 jstart = halo + 1
 iend   = i_target - halo
 jend   = j_target - halo

 allocate(lsoil_data(lsoil_target))
 do i = 1, lsoil_target
   lsoil_data(i) = float(i)
 enddo

 allocate(x_data(i_target_out))
 do i = 1, i_target_out
   x_data(i) = float(i)
 enddo

 allocate(y_data(j_target_out))
 do i = 1, j_target_out
   y_data(i) = float(i)
 enddo

 if (convert_nst) then
   print*,'- WRITE FV3 SURFACE AND NST DATA TO NETCDF FILE'
 else
   print*,'- WRITE FV3 SURFACE DATA TO NETCDF FILE'
 endif

 if (localpet == 0) then
   allocate(data_one_tile(i_target,j_target))
   allocate(data_one_tile_3d(i_target,j_target,lsoil_target))
   allocate(idata_one_tile(i_target,j_target))
   allocate(dum2d(i_target_out,j_target_out))
   allocate(dum3d(i_target_out,j_target_out,lsoil_target))
 else
   allocate(data_one_tile(0,0))
   allocate(data_one_tile_3d(0,0,0))
   allocate(idata_one_tile(0,0))
   allocate(dum2d(0,0))
   allocate(dum3d(0,0,0))
 endif

 TILE_LOOP : do tile = 1, num_tiles_target_grid

   LOCAL_PET : if (localpet == 0) then

!-- SGregory mod... sfc filename     WRITE(OUTFILE, '(A, I1, A)'), 'out.sfc.tile', tile, '.nc'
     WRITE(OUTFILE, '(A, I1, A)'), 'sfc_data.tile', tile, '.nc'
!--- open the file
     error = nf90_create(outfile, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                         ncid, initialsize=inital, chunksize=fsize)
     call netcdf_err(error, 'CREATING FILE='//trim(outfile) )

!--- define dimensions
     error = nf90_def_dim(ncid, 'xaxis_1', i_target_out, dim_x)
     call netcdf_err(error, 'DEFINING XAXIS DIMENSION' )
     error = nf90_def_dim(ncid, 'yaxis_1', j_target_out, dim_y)
     call netcdf_err(error, 'DEFINING YAXIS DIMENSION' )
     error = nf90_def_dim(ncid, 'zaxis_1', lsoil_target, dim_lsoil)
     call netcdf_err(error, 'DEFINING ZAXIS DIMENSION' )
     error = nf90_def_dim(ncid, 'Time', 1, dim_time)
     call netcdf_err(error, 'DEFINING TIME DIMENSION' )

 !--- define fields
     error = nf90_def_var(ncid, 'xaxis_1', NF90_FLOAT, (/dim_x/), id_x)
     call netcdf_err(error, 'DEFINING XAXIS_1 FIELD' )
     error = nf90_put_att(ncid, id_x, "long_name", "xaxis_1")
     call netcdf_err(error, 'DEFINING XAXIS_1 LONG NAME' )
     error = nf90_put_att(ncid, id_x, "units", "none")
     call netcdf_err(error, 'DEFINING XAXIS_1 UNITS' )
     error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
     call netcdf_err(error, 'WRITING XAXIS_1 FIELD' )

     error = nf90_def_var(ncid, 'yaxis_1', NF90_FLOAT, (/dim_y/), id_y)
     call netcdf_err(error, 'DEFINING YAXIS_1 FIELD' )
     error = nf90_put_att(ncid, id_y, "long_name", "yaxis_1")
     call netcdf_err(error, 'DEFINING YAXIS_1 LONG NAME' )
     error = nf90_put_att(ncid, id_y, "units", "none")
     call netcdf_err(error, 'DEFINING YAXIS_1 UNITS' )
     error = nf90_put_att(ncid, id_y, "cartesian_axis", "Y")
     call netcdf_err(error, 'WRITING YAXIS_1 FIELD' )

     error = nf90_def_var(ncid, 'zaxis_1', NF90_FLOAT, (/dim_lsoil/), id_lsoil)
     call netcdf_err(error, 'DEFINING ZAXIS_1 FIELD' )
     error = nf90_put_att(ncid, id_lsoil, "long_name", "zaxis_1")
     call netcdf_err(error, 'DEFINING ZAXIS_1 LONG NAME' )
     error = nf90_put_att(ncid, id_lsoil, "units", "none")
     call netcdf_err(error, 'DEFINING ZAXIS_1 UNITS' )
     error = nf90_put_att(ncid, id_lsoil, "cartesian_axis", "Z")
     call netcdf_err(error, 'WRITING ZAXIS_1 FIELD' )

     error = nf90_def_var(ncid, 'Time', NF90_FLOAT, dim_time, id_time)
     call netcdf_err(error, 'DEFINING TIME FIELD' )
     error = nf90_put_att(ncid, id_time, "long_name", "Time")
     call netcdf_err(error, 'DEFINING TIME LONG NAME' )
     error = nf90_put_att(ncid, id_time, "units", "time level")
     call netcdf_err(error, 'DEFINING TIME UNITS' )
     error = nf90_put_att(ncid, id_time, "cartesian_axis", "T")
     call netcdf_err(error, 'WRITING TIME FIELD' )

     error = nf90_def_var(ncid, 'slmsk', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_slmsk)
     call netcdf_err(error, 'DEFINING SLMSK' )
     error = nf90_put_att(ncid, id_slmsk, "long_name", "slmsk")
     call netcdf_err(error, 'DEFINING SLMSK LONG NAME' )
     error = nf90_put_att(ncid, id_slmsk, "units", "none")
     call netcdf_err(error, 'DEFINING SLMSK UNITS' )

     error = nf90_def_var(ncid, 'tsea', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_tsea)
     call netcdf_err(error, 'DEFINING TSEA' )
     error = nf90_put_att(ncid, id_tsea, "long_name", "tsea")
     call netcdf_err(error, 'DEFINING TSEA LONG NAME' )
     error = nf90_put_att(ncid, id_tsea, "units", "none")
     call netcdf_err(error, 'DEFINING TSEA UNITS' )

     error = nf90_def_var(ncid, 'sheleg', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_sheleg)
     call netcdf_err(error, 'DEFINING SHELEG' )
     error = nf90_put_att(ncid, id_sheleg, "long_name", "sheleg")
     call netcdf_err(error, 'DEFINING SHELEG LONG NAME' )
     error = nf90_put_att(ncid, id_sheleg, "units", "none")
     call netcdf_err(error, 'DEFINING SHELEG UNITS' )

     error = nf90_def_var(ncid, 'tg3', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_tg3)
     call netcdf_err(error, 'DEFINING TG3' )
     error = nf90_put_att(ncid, id_tg3, "long_name", "tg3")
     call netcdf_err(error, 'DEFINING TG3 LONG NAME' )
     error = nf90_put_att(ncid, id_tg3, "units", "none")
     call netcdf_err(error, 'DEFINING TG3 UNITS' )

     error = nf90_def_var(ncid, 'zorl', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_zorl)
     call netcdf_err(error, 'DEFINING ZORL' )
     error = nf90_put_att(ncid, id_zorl, "long_name", "zorl")
     call netcdf_err(error, 'DEFINING ZORL LONG NAME' )
     error = nf90_put_att(ncid, id_zorl, "units", "none")
     call netcdf_err(error, 'DEFINING ZORL UNITS' )

     error = nf90_def_var(ncid, 'alvsf', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_alvsf)
     call netcdf_err(error, 'DEFINING ALVSF' )
     error = nf90_put_att(ncid, id_alvsf, "long_name", "alvsf")
     call netcdf_err(error, 'DEFINING ALVSF LONG NAME' )
     error = nf90_put_att(ncid, id_alvsf, "units", "none")
     call netcdf_err(error, 'DEFINING ALVSF UNITS' )

     error = nf90_def_var(ncid, 'alvwf', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_alvwf)
     call netcdf_err(error, 'DEFINING ALVWF' )
     error = nf90_put_att(ncid, id_alvwf, "long_name", "alvwf")
     call netcdf_err(error, 'DEFINING ALVWF LONG NAME' )
     error = nf90_put_att(ncid, id_alvwf, "units", "none")
     call netcdf_err(error, 'DEFINING ALVWF UNITS' )

     error = nf90_def_var(ncid, 'alnsf', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_alnsf)
     call netcdf_err(error, 'DEFINING ALNSF' )
     error = nf90_put_att(ncid, id_alnsf, "long_name", "alnsf")
     call netcdf_err(error, 'DEFINING ALNSF LONG NAME' )
     error = nf90_put_att(ncid, id_alnsf, "units", "none")
     call netcdf_err(error, 'DEFINING ALNSF UNITS' )

     error = nf90_def_var(ncid, 'alnwf', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_alnwf)
     call netcdf_err(error, 'DEFINING ALNWF' )
     error = nf90_put_att(ncid, id_alnwf, "long_name", "alnwf")
     call netcdf_err(error, 'DEFINING ALNWF LONG NAME' )
     error = nf90_put_att(ncid, id_alnwf, "units", "none")
     call netcdf_err(error, 'DEFINING ALNWF UNITS' )

     error = nf90_def_var(ncid, 'facsf', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_facsf)
     call netcdf_err(error, 'DEFINING FACSF' )
     error = nf90_put_att(ncid, id_facsf, "long_name", "facsf")
     call netcdf_err(error, 'DEFINING FACSF LONG NAME' )
     error = nf90_put_att(ncid, id_facsf, "units", "none")
     call netcdf_err(error, 'DEFINING FACSF UNITS' )

     error = nf90_def_var(ncid, 'facwf', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_facwf)
     call netcdf_err(error, 'DEFINING FACWF' )
     error = nf90_put_att(ncid, id_facwf, "long_name", "facwf")
     call netcdf_err(error, 'DEFINING FACWF LONG NAME' )
     error = nf90_put_att(ncid, id_facwf, "units", "none")
     call netcdf_err(error, 'DEFINING FACWF UNITS' )

     error = nf90_def_var(ncid, 'vfrac', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_vfrac)
     call netcdf_err(error, 'DEFINING VFRAC' )
     error = nf90_put_att(ncid, id_vfrac, "long_name", "vfrac")
     call netcdf_err(error, 'DEFINING VFRAC LONG NAME' )
     error = nf90_put_att(ncid, id_vfrac, "units", "none")
     call netcdf_err(error, 'DEFINING VFRAC UNITS' )

     error = nf90_def_var(ncid, 'canopy', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_canopy)
     call netcdf_err(error, 'DEFINING CANOPY' )
     error = nf90_put_att(ncid, id_canopy, "long_name", "canopy")
     call netcdf_err(error, 'DEFINING CANOPY LONG NAME' )
     error = nf90_put_att(ncid, id_canopy, "units", "none")
     call netcdf_err(error, 'DEFINING CANOPY UNITS' )

     error = nf90_def_var(ncid, 'f10m', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_f10m)
     call netcdf_err(error, 'DEFINING F10M' )
     error = nf90_put_att(ncid, id_f10m, "long_name", "f10m")
     call netcdf_err(error, 'DEFINING F10M LONG NAME' )
     error = nf90_put_att(ncid, id_f10m, "units", "none")
     call netcdf_err(error, 'DEFINING F10M UNITS' )

     error = nf90_def_var(ncid, 't2m', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_t2m)
     call netcdf_err(error, 'DEFINING T2M' )
     error = nf90_put_att(ncid, id_t2m, "long_name", "t2m")
     call netcdf_err(error, 'DEFINING T2M LONG NAME' )
     error = nf90_put_att(ncid, id_t2m, "units", "none")
     call netcdf_err(error, 'DEFINING T2M UNITS' )

     error = nf90_def_var(ncid, 'q2m', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_q2m)
     call netcdf_err(error, 'DEFINING Q2M' )
     error = nf90_put_att(ncid, id_q2m, "long_name", "q2m")
     call netcdf_err(error, 'DEFINING Q2M LONG NAME' )
     error = nf90_put_att(ncid, id_q2m, "units", "none")
     call netcdf_err(error, 'DEFINING Q2M UNITS' )

     error = nf90_def_var(ncid, 'vtype', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_vtype)
     call netcdf_err(error, 'DEFINING VTYPE' )
     error = nf90_put_att(ncid, id_vtype, "long_name", "vtype")
     call netcdf_err(error, 'DEFINING VTYPE LONG NAME' )
     error = nf90_put_att(ncid, id_vtype, "units", "none")
     call netcdf_err(error, 'DEFINING VTYPE UNITS' )

     error = nf90_def_var(ncid, 'stype', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_stype)
     call netcdf_err(error, 'DEFINING STYPE' )
     error = nf90_put_att(ncid, id_stype, "long_name", "stype")
     call netcdf_err(error, 'DEFINING STYPE LONG NAME' )
     error = nf90_put_att(ncid, id_stype, "units", "none")
     call netcdf_err(error, 'DEFINING STYPE UNITS' )

     error = nf90_def_var(ncid, 'uustar', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_uustar)
     call netcdf_err(error, 'DEFINING UUSTAR' )
     error = nf90_put_att(ncid, id_uustar, "long_name", "uustar")
     call netcdf_err(error, 'DEFINING UUSTAR LONG NAME' )
     error = nf90_put_att(ncid, id_uustar, "units", "none")
     call netcdf_err(error, 'DEFINING UUSTAR UNITS' )

     error = nf90_def_var(ncid, 'ffmm', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_ffmm)
     call netcdf_err(error, 'DEFINING FFMM' )
     error = nf90_put_att(ncid, id_ffmm, "long_name", "ffmm")
     call netcdf_err(error, 'DEFINING FFMM LONG NAME' )
     error = nf90_put_att(ncid, id_ffmm, "units", "none")
     call netcdf_err(error, 'DEFINING FFMM UNITS' )

     error = nf90_def_var(ncid, 'ffhh', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_ffhh)
     call netcdf_err(error, 'DEFINING FFHH' )
     error = nf90_put_att(ncid, id_ffhh, "long_name", "ffhh")
     call netcdf_err(error, 'DEFINING FFHH LONG NAME' )
     error = nf90_put_att(ncid, id_ffhh, "units", "none")
     call netcdf_err(error, 'DEFINING FFHH UNITS' )

     error = nf90_def_var(ncid, 'hice', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_hice)
     call netcdf_err(error, 'DEFINING HICE' )
     error = nf90_put_att(ncid, id_hice, "long_name", "hice")
     call netcdf_err(error, 'DEFINING HICE LONG NAME' )
     error = nf90_put_att(ncid, id_hice, "units", "none")
     call netcdf_err(error, 'DEFINING HICE UNITS' )

     error = nf90_def_var(ncid, 'fice', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_fice)
     call netcdf_err(error, 'DEFINING FICE' )
     error = nf90_put_att(ncid, id_fice, "long_name", "fice")
     call netcdf_err(error, 'DEFINING FICE LONG NAME' )
     error = nf90_put_att(ncid, id_fice, "units", "none")
     call netcdf_err(error, 'DEFINING FICE UNITS' )

     error = nf90_def_var(ncid, 'tisfc', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_tisfc)
     call netcdf_err(error, 'DEFINING TISFC' )
     error = nf90_put_att(ncid, id_tisfc, "long_name", "tisfc")
     call netcdf_err(error, 'DEFINING TISFC LONG NAME' )
     error = nf90_put_att(ncid, id_tisfc, "units", "none")
     call netcdf_err(error, 'DEFINING TISFC UNITS' )

     error = nf90_def_var(ncid, 'tprcp', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_tprcp)
     call netcdf_err(error, 'DEFINING TPRCP' )
     error = nf90_put_att(ncid, id_tprcp, "long_name", "tprcp")
     call netcdf_err(error, 'DEFINING TPRCP LONG NAME' )
     error = nf90_put_att(ncid, id_tprcp, "units", "none")
     call netcdf_err(error, 'DEFINING TPRCP UNITS' )

     error = nf90_def_var(ncid, 'srflag', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_srflag)
     call netcdf_err(error, 'DEFINING SRFLAG' )
     error = nf90_put_att(ncid, id_srflag, "long_name", "srflag")
     call netcdf_err(error, 'DEFINING SRFLAG LONG NAME' )
     error = nf90_put_att(ncid, id_srflag, "units", "none")
     call netcdf_err(error, 'DEFINING SRFLAG UNITS' )

     error = nf90_def_var(ncid, 'snwdph', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_snwdph)
     call netcdf_err(error, 'DEFINING SNWDPH' )
     error = nf90_put_att(ncid, id_snwdph, "long_name", "snwdph")
     call netcdf_err(error, 'DEFINING SNWDPH LONG NAME' )
     error = nf90_put_att(ncid, id_snwdph, "units", "none")
     call netcdf_err(error, 'DEFINING SNWDPH UNITS' )

     error = nf90_def_var(ncid, 'shdmin', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_shdmin)
     call netcdf_err(error, 'DEFINING SHDMIN' )
     error = nf90_put_att(ncid, id_shdmin, "long_name", "shdmin")
     call netcdf_err(error, 'DEFINING SHDMIN LONG NAME' )
     error = nf90_put_att(ncid, id_shdmin, "units", "none")
     call netcdf_err(error, 'DEFINING SHDMIN UNITS' )

     error = nf90_def_var(ncid, 'shdmax', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_shdmax)
     call netcdf_err(error, 'DEFINING SHDMAX' )
     error = nf90_put_att(ncid, id_shdmax, "long_name", "shdmax")
     call netcdf_err(error, 'DEFINING SHDMAX LONG NAME' )
     error = nf90_put_att(ncid, id_shdmax, "units", "none")
     call netcdf_err(error, 'DEFINING SHDMAX UNITS' )

     error = nf90_def_var(ncid, 'slope', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_slope)
     call netcdf_err(error, 'DEFINING SLOPE' )
     error = nf90_put_att(ncid, id_slope, "long_name", "slope")
     call netcdf_err(error, 'DEFINING SLOPE LONG NAME' )
     error = nf90_put_att(ncid, id_slope, "units", "none")
     call netcdf_err(error, 'DEFINING SLOPE UNITS' )

     error = nf90_def_var(ncid, 'snoalb', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_snoalb)
     call netcdf_err(error, 'DEFINING SNOALB' )
     error = nf90_put_att(ncid, id_snoalb, "long_name", "snoalb")
     call netcdf_err(error, 'DEFINING SNOALB LONG NAME' )
     error = nf90_put_att(ncid, id_snoalb, "units", "none")
     call netcdf_err(error, 'DEFINING SNOALB UNITS' )

     error = nf90_def_var(ncid, 'stc', NF90_DOUBLE, (/dim_x,dim_y,dim_lsoil,dim_time/), id_stc)
     call netcdf_err(error, 'DEFINING STC' )
     error = nf90_put_att(ncid, id_stc, "long_name", "stc")
     call netcdf_err(error, 'DEFINING STC LONG NAME' )
     error = nf90_put_att(ncid, id_stc, "units", "none")
     call netcdf_err(error, 'DEFINING STC UNITS' )

     error = nf90_def_var(ncid, 'smc', NF90_DOUBLE, (/dim_x,dim_y,dim_lsoil,dim_time/), id_smc)
     call netcdf_err(error, 'DEFINING SMC' )
     error = nf90_put_att(ncid, id_smc, "long_name", "smc")
     call netcdf_err(error, 'DEFINING SMC LONG NAME' )
     error = nf90_put_att(ncid, id_smc, "units", "none")
     call netcdf_err(error, 'DEFINING SMC UNITS' )

     error = nf90_def_var(ncid, 'slc', NF90_DOUBLE, (/dim_x,dim_y,dim_lsoil,dim_time/), id_slc)
     call netcdf_err(error, 'DEFINING SLC' )
     error = nf90_put_att(ncid, id_slc, "long_name", "slc")
     call netcdf_err(error, 'DEFINING SLC LONG NAME' )
     error = nf90_put_att(ncid, id_slc, "units", "none")
     call netcdf_err(error, 'DEFINING SLC UNITS' )

     if (convert_nst) then

       error = nf90_def_var(ncid, 'tref', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_tref)
       call netcdf_err(error, 'DEFINING TREF' )
       error = nf90_put_att(ncid, id_tref, "long_name", "tref")
       call netcdf_err(error, 'DEFINING TREF LONG NAME' )
       error = nf90_put_att(ncid, id_tref, "units", "none")
       call netcdf_err(error, 'DEFINING TREF UNITS' )

       error = nf90_def_var(ncid, 'z_c', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_z_c)
       call netcdf_err(error, 'DEFINING Z_C' )
       error = nf90_put_att(ncid, id_z_c, "long_name", "z_c")
       call netcdf_err(error, 'DEFINING Z_C LONG NAME' )
       error = nf90_put_att(ncid, id_z_c, "units", "none")
       call netcdf_err(error, 'DEFINING Z_C UNITS' )

       error = nf90_def_var(ncid, 'c_0', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_c_0)
       call netcdf_err(error, 'DEFINING C_0' )
       error = nf90_put_att(ncid, id_c_0, "long_name", "c_0")
       call netcdf_err(error, 'DEFINING C_0 LONG NAME' )
       error = nf90_put_att(ncid, id_c_0, "units", "none")
       call netcdf_err(error, 'DEFINING C_0 UNITS' )

       error = nf90_def_var(ncid, 'c_d', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_c_d)
       call netcdf_err(error, 'DEFINING C_D' )
       error = nf90_put_att(ncid, id_c_d, "long_name", "c_d")
       call netcdf_err(error, 'DEFINING C_D LONG NAME' )
       error = nf90_put_att(ncid, id_c_d, "units", "none")
       call netcdf_err(error, 'DEFINING C_D UNITS' )

       error = nf90_def_var(ncid, 'w_0', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_w_0)
       call netcdf_err(error, 'DEFINING W_0' )
       error = nf90_put_att(ncid, id_w_0, "long_name", "w_0")
       call netcdf_err(error, 'DEFINING W_0 LONG NAME' )
       error = nf90_put_att(ncid, id_w_0, "units", "none")
       call netcdf_err(error, 'DEFINING W_0 UNITS' )

       error = nf90_def_var(ncid, 'w_d', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_w_d)
       call netcdf_err(error, 'DEFINING W_D' )
       error = nf90_put_att(ncid, id_w_d, "long_name", "w_d")
       call netcdf_err(error, 'DEFINING W_D LONG NAME' )
       error = nf90_put_att(ncid, id_w_d, "units", "none")
       call netcdf_err(error, 'DEFINING W_D UNITS' )

       error = nf90_def_var(ncid, 'xt', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_xt)
       call netcdf_err(error, 'DEFINING XT' )
       error = nf90_put_att(ncid, id_xt, "long_name", "xt")
       call netcdf_err(error, 'DEFINING XT LONG NAME' )
       error = nf90_put_att(ncid, id_xt, "units", "none")
       call netcdf_err(error, 'DEFINING XT UNITS' )

       error = nf90_def_var(ncid, 'xs', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_xs)
       call netcdf_err(error, 'DEFINING XS' )
       error = nf90_put_att(ncid, id_xs, "long_name", "xs")
       call netcdf_err(error, 'DEFINING XS LONG NAME' )
       error = nf90_put_att(ncid, id_xs, "units", "none")
       call netcdf_err(error, 'DEFINING XS UNITS' )

       error = nf90_def_var(ncid, 'xu', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_xu)
       call netcdf_err(error, 'DEFINING XU' )
       error = nf90_put_att(ncid, id_xu, "long_name", "xu")
       call netcdf_err(error, 'DEFINING XU LONG NAME' )
       error = nf90_put_att(ncid, id_xu, "units", "none")
       call netcdf_err(error, 'DEFINING XU UNITS' )

       error = nf90_def_var(ncid, 'xv', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_xv)
       call netcdf_err(error, 'DEFINING XV' )
       error = nf90_put_att(ncid, id_xv, "long_name", "xv")
       call netcdf_err(error, 'DEFINING XV LONG NAME' )
       error = nf90_put_att(ncid, id_xv, "units", "none")
       call netcdf_err(error, 'DEFINING XV UNITS' )

       error = nf90_def_var(ncid, 'xz', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_xz)
       call netcdf_err(error, 'DEFINING XZ' )
       error = nf90_put_att(ncid, id_xz, "long_name", "xz")
       call netcdf_err(error, 'DEFINING XZ LONG NAME' )
       error = nf90_put_att(ncid, id_xz, "units", "none")
       call netcdf_err(error, 'DEFINING XZ UNITS' )

       error = nf90_def_var(ncid, 'zm', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_zm)
       call netcdf_err(error, 'DEFINING ZM' )
       error = nf90_put_att(ncid, id_zm, "long_name", "zm")
       call netcdf_err(error, 'DEFINING ZM LONG NAME' )
       error = nf90_put_att(ncid, id_zm, "units", "none")
       call netcdf_err(error, 'DEFINING ZM UNITS' )

       error = nf90_def_var(ncid, 'xtts', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_xtts)
       call netcdf_err(error, 'DEFINING XTTS' )
       error = nf90_put_att(ncid, id_xtts, "long_name", "xtts")
       call netcdf_err(error, 'DEFINING XTTS LONG NAME' )
       error = nf90_put_att(ncid, id_xtts, "units", "none")
       call netcdf_err(error, 'DEFINING XTTS UNITS' )

       error = nf90_def_var(ncid, 'xzts', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_xzts)
       call netcdf_err(error, 'DEFINING XZTS' )
       error = nf90_put_att(ncid, id_xzts, "long_name", "xzts")
       call netcdf_err(error, 'DEFINING XZTS LONG NAME' )
       error = nf90_put_att(ncid, id_xzts, "units", "none")
       call netcdf_err(error, 'DEFINING XZTS UNITS' )

       error = nf90_def_var(ncid, 'd_conv', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_d_conv)
       call netcdf_err(error, 'DEFINING D_CONV' )
       error = nf90_put_att(ncid, id_d_conv, "long_name", "d_conv")
       call netcdf_err(error, 'DEFINING D_CONV LONG NAME' )
       error = nf90_put_att(ncid, id_d_conv, "units", "none")
       call netcdf_err(error, 'DEFINING D_CONV UNITS' )

       error = nf90_def_var(ncid, 'ifd', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_ifd)
       call netcdf_err(error, 'DEFINING IFD' )
       error = nf90_put_att(ncid, id_ifd, "long_name", "ifd")
       call netcdf_err(error, 'DEFINING IFD LONG NAME' )
       error = nf90_put_att(ncid, id_ifd, "units", "none")
       call netcdf_err(error, 'DEFINING IFD UNITS' )

       error = nf90_def_var(ncid, 'dt_cool', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_dt_cool)
       call netcdf_err(error, 'DEFINING DT_COOL' )
       error = nf90_put_att(ncid, id_dt_cool, "long_name", "dt_cool")
       call netcdf_err(error, 'DEFINING DT_COOL LONG NAME' )
       error = nf90_put_att(ncid, id_dt_cool, "units", "none")
       call netcdf_err(error, 'DEFINING DT_COOL UNITS' )

       error = nf90_def_var(ncid, 'qrain', NF90_DOUBLE, (/dim_x,dim_y,dim_time/), id_qrain)
       call netcdf_err(error, 'DEFINING QRAIN' )
       error = nf90_put_att(ncid, id_qrain, "long_name", "qrain")
       call netcdf_err(error, 'DEFINING QRAIN LONG NAME' )
       error = nf90_put_att(ncid, id_qrain, "units", "none")
       call netcdf_err(error, 'DEFINING QRAIN UNITS' )

     endif  ! nsst records

     error = nf90_enddef(ncid, header_buffer_val,4,0,4)
     call netcdf_err(error, 'DEFINING HEADER' )

   endif LOCAL_PET ! is localpet 0?

   if (localpet == 0) then
     error = nf90_put_var( ncid, id_lsoil, lsoil_data)
     call netcdf_err(error, 'WRITING ZAXIS RECORD' )
     error = nf90_put_var( ncid, id_x, x_data)
     call netcdf_err(error, 'WRITING XAXIS RECORD' )
     error = nf90_put_var( ncid, id_y, y_data)
     call netcdf_err(error, 'WRITING YAXIS RECORD' )
     times = 1.0
     error = nf90_put_var( ncid, id_time, times)
     call netcdf_err(error, 'WRITING TIME RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SNOW LIQ EQUIV FOR TILE: ", tile
   call ESMF_FieldGather(snow_liq_equiv_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_sheleg, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING SNOW LIQ EQUIV RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SNOW DEPTH FOR TILE: ", tile
   call ESMF_FieldGather(snow_depth_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_snwdph, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING SNWDPH RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SLOPE TYPE FOR TILE: ", tile
   call ESMF_FieldGather(slope_type_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_slope, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING SLOPE RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID Z0 FOR TILE: ", tile
   call ESMF_FieldGather(z0_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_zorl, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING Z0 RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID MAX SNOW ALBEDO FOR TILE: ", tile
   call ESMF_FieldGather(mxsno_albedo_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_snoalb, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING MAX SNOW ALBEDO RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SOIL TYPE FOR TILE: ", tile
   call ESMF_FieldGather(soil_type_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_stype, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING SOIL TYPE RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID VEGETATION TYPE FOR TILE: ", tile
   call ESMF_FieldGather(veg_type_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_vtype, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING VEGETATION TYPE RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID VEGETATION GREENNESS FOR TILE: ", tile
   call ESMF_FieldGather(veg_greenness_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_vfrac, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING VEGETATION GREENNESS RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SUBSTRATE TEMPERATURE FOR TILE: ", tile
   call ESMF_FieldGather(substrate_temp_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_tg3, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING SUBSTRATE TEMPERATURE RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID FACSF FOR TILE: ", tile
   call ESMF_FieldGather(facsf_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_facsf, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING FACSF RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID FACWF FOR TILE: ", tile
   call ESMF_FieldGather(facwf_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_facwf, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING FACWF RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID ALNSF FOR TILE: ", tile
   call ESMF_FieldGather(alnsf_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_alnsf, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING ALNSF RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID ALNWF FOR TILE: ", tile
   call ESMF_FieldGather(alnwf_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_alnwf, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING ALNWF RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID ALVSF FOR TILE: ", tile
   call ESMF_FieldGather(alvsf_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_alvsf, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING ALVSF RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID ALVWF FOR TILE: ", tile
   call ESMF_FieldGather(alvwf_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_alvwf, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING ALVWF RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID MAX VEGETATION GREENNESS FOR TILE: ", tile
   call ESMF_FieldGather(max_veg_greenness_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_shdmax, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING MAX VEGETATION GREENNESS RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID MIN VEGETATION GREENNESS FOR TILE: ", tile
   call ESMF_FieldGather(min_veg_greenness_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_shdmin, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING MIN VEGETATION GREENNESS RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID T2M FOR TILE: ", tile
   call ESMF_FieldGather(t2m_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_t2m, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING T2M RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID Q2M FOR TILE: ", tile
   call ESMF_FieldGather(q2m_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_q2m, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING Q2M RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID TPRCP FOR TILE: ", tile
   call ESMF_FieldGather(tprcp_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_tprcp, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING TPRCP RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID F10M FOR TILE: ", tile
   call ESMF_FieldGather(f10m_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_f10m, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING F10M RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID FFMM FOR TILE: ", tile
   call ESMF_FieldGather(ffmm_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_ffmm, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING FFMM RECORD' )
     dum2d = 0.0
     error = nf90_put_var( ncid, id_ffhh, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING FFHH RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID USTAR FOR TILE: ", tile
   call ESMF_FieldGather(ustar_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_uustar, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING USTAR RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SRFLAG FOR TILE: ", tile
   call ESMF_FieldGather(srflag_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_srflag, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING SRFLAG RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SEA ICE FRACTION FOR TILE: ", tile
   call ESMF_FieldGather(seaice_fract_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_fice, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING FICE RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SEA ICE DEPTH FOR TILE: ", tile
   call ESMF_FieldGather(seaice_depth_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_hice, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING HICE RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SEA ICE SKIN TEMP FOR TILE: ", tile
   call ESMF_FieldGather(seaice_skin_temp_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_tisfc, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING TISFC RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SKIN TEMP FOR TILE: ", tile
   call ESMF_FieldGather(skin_temp_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_tsea, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING TSEA RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID LANDMASK FOR TILE: ", tile
   call ESMF_FieldGather(landmask_target_grid, idata_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = float(idata_one_tile(istart:iend, jstart:jend))
     error = nf90_put_var( ncid, id_slmsk, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING LANDMASK RECORD' )
   endif

   print*,"- CALL FieldGather FOR TARGET GRID CANOPY MOISTURE CONTENT FOR TILE: ", tile
   call ESMF_FieldGather(canopy_mc_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
     error = nf90_put_var( ncid, id_canopy, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
     call netcdf_err(error, 'WRITING CANOPY MC RECORD' )
   endif

! soil temperature 

   print*,"- CALL FieldGather FOR TARGET GRID SOIL TEMPERATURE FOR TILE: ", tile
   call ESMF_FieldGather(soil_temp_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum3d(:,:,:) = data_one_tile_3d(istart:iend, jstart:jend,:)
     error = nf90_put_var( ncid, id_stc, dum3d, start=(/1,1,1,1/), count=(/i_target_out,j_target_out,lsoil_target,1/))
     call netcdf_err(error, 'WRITING SOIL TEMP RECORD' )
   endif

! soil moisture (total)

   print*,"- CALL FieldGather FOR TARGET GRID TOTAL SOIL MOISTURE FOR TILE: ", tile
   call ESMF_FieldGather(soilm_tot_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum3d(:,:,:) = data_one_tile_3d(istart:iend, jstart:jend,:)
     error = nf90_put_var( ncid, id_smc, dum3d, start=(/1,1,1,1/), count=(/i_target_out,j_target_out,lsoil_target,1/))
     call netcdf_err(error, 'WRITING TOTAL SOIL MOISTURE RECORD' )
   endif

! soil moisture (liquid)

   print*,"- CALL FieldGather FOR TARGET GRID LIQUID SOIL MOISTURE FOR TILE: ", tile
   call ESMF_FieldGather(soilm_liq_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call mpi_abort

   if (localpet == 0) then
     dum3d(:,:,:) = data_one_tile_3d(istart:iend, jstart:jend,:)
     error = nf90_put_var( ncid, id_slc, dum3d, start=(/1,1,1,1/), count=(/i_target_out,j_target_out,lsoil_target,1/))
     call netcdf_err(error, 'WRITING LIQUID SOIL MOISTURE RECORD' )
   endif

   if (convert_nst) then

     print*,"- CALL FieldGather FOR TARGET C_D FOR TILE: ", tile
     call ESMF_FieldGather(c_d_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_c_d, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING C_D RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET C_0 FOR TILE: ", tile
     call ESMF_FieldGather(c_0_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_c_0, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING C_0 RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET D_CONV FOR TILE: ", tile
     call ESMF_FieldGather(d_conv_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_d_conv, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING D_CONV RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET DT_COOL FOR TILE: ", tile
     call ESMF_FieldGather(dt_cool_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_dt_cool, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING DT_COOL RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET IFD FOR TILE: ", tile
     call ESMF_FieldGather(ifd_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_ifd, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING IFD RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET QRAIN FOR TILE: ", tile
     call ESMF_FieldGather(qrain_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_qrain, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING QRAIN RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET TREF FOR TILE: ", tile
     call ESMF_FieldGather(tref_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_tref, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING TREF RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET W_D FOR TILE: ", tile
     call ESMF_FieldGather(w_d_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_w_d, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING W_D RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET W_0 FOR TILE: ", tile
     call ESMF_FieldGather(w_0_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_w_0, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING W_0 RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET XS FOR TILE: ", tile
     call ESMF_FieldGather(xs_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_xs, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING XS RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET XT FOR TILE: ", tile
     call ESMF_FieldGather(xt_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_xt, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING XT RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET XU FOR TILE: ", tile
     call ESMF_FieldGather(xu_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_xu, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING XU RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET XV FOR TILE: ", tile
     call ESMF_FieldGather(xv_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_xv, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING XV RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET XZ FOR TILE: ", tile
     call ESMF_FieldGather(xz_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_xz, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING XZ RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET XTTS FOR TILE: ", tile
     call ESMF_FieldGather(xtts_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_xtts, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING XTTS RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET XZTS FOR TILE: ", tile
     call ESMF_FieldGather(xzts_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_xzts, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING XZTS RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET Z_C FOR TILE: ", tile
     call ESMF_FieldGather(z_c_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_z_c, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING Z_C RECORD' )
     endif

     print*,"- CALL FieldGather FOR TARGET ZM FOR TILE: ", tile
     call ESMF_FieldGather(zm_target_grid, data_one_tile, rootPet=0, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
        call mpi_abort

     if (localpet == 0) then
       dum2d(:,:) = data_one_tile(istart:iend, jstart:jend)
       error = nf90_put_var( ncid, id_zm, dum2d, start=(/1,1,1/), count=(/i_target_out,j_target_out,1/))
       call netcdf_err(error, 'WRITING ZM RECORD' )
     endif

   endif ! convert nst

!-------------------------------------------------------------------------------
! close file
!-------------------------------------------------------------------------------

   error = nf90_close(ncid)

 enddo TILE_LOOP

 deallocate(lsoil_data, x_data, y_data)
 deallocate(data_one_tile, data_one_tile_3d, idata_one_tile, dum2d, dum3d)

 return

 end subroutine write_fv3_sfc_data_netcdf
