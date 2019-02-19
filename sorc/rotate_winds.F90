program rotate_winds
use NETCDF
implicit none


integer  :: i,j,k,nx,ny,CRES,km

integer :: netid,stat,varid,tile
real inner_prod
character*1 :: tilef
real,allocatable,dimension(:) :: p1,p2,p3,ex,ey,e1
real,allocatable,dimension(:,:)  :: lat,lon,lat2,lon2
real,allocatable,dimension(:,:,:)  :: grid3,ec1,ec2,vlat,vlon
real,allocatable,dimension(:,:,:)  ::
ud,vd,ua,va,udnew,vdnew,u_s,v_s,u_w,v_w

! allocate arrays
km=1
CRES=128
nx=CRES
ny=CRES
! fill lat-lon definition
allocate(lon2(2*nx+1,2*ny+1))
allocate(lat2(2*nx+1,2*ny+1))
! b-grid grid definition
allocate(lat(nx+1,ny+1))
allocate(lon(nx+1,ny+1))

! wind fields
allocate(ud(nx  ,ny+1,km))
allocate(vd(nx+1,ny ,km ))
allocate(udnew(nx  ,ny+1,km))
allocate(vdnew(nx+1,ny ,km ))
allocate(ua(nx,ny,km))
allocate(va(nx,ny,km))
allocate(u_s(nx  ,ny+1,km))
allocate(v_s(nx  ,ny+1,km))
allocate(u_w(nx+1,ny ,km ))
allocate(v_w(nx+1,ny  ,km))
! vector stuff
allocate(p1(2))
allocate(p2(2))
allocate(p3(2))
allocate(e1(3))
allocate(ex(3))
allocate(ey(3))

! open grid file
do tile=1,6
write(tilef,FMT='(I1)') tile
print*,'tilef'
stat=nf90_open('/scratch3/BMC/gsienkf/Philip.Pegion/FV3/experiments/C96/stoch_master2/2014080100/mem01/INPUT/C96_grid.tile'//tilef//'.nc',NF90_NOWRITE,netid)
stat=nf90_inq_varid(netid,'x',varid)
stat=nf90_get_var(netid,varid,lon2,count=(/2*nx+1,2*ny+1/))
stat=nf90_inq_varid(netid,'y',varid)
stat=nf90_get_var(netid,varid,lat2,count=(/2*nx+1,2*ny+1/))
! convert latitude and longitude to radians
lon2=lon2*3.1415928/180.0
lat2=lat2*3.1415928/180.0
print*,'lat=',maxval(lat2),minval(lat2)
print*,'lon=',maxval(lon2),minval(lon2)
stat=nf90_close(netid)

! open restart_file
stat=nf90_open('/scratch3/BMC/gsienkf/Philip.Pegion/FV3/experiments/C96/stoch_master2/2014080100/mem01/RESTART/fv_core.res.tile'//tilef//'.nc',NF90_NOWRITE,netid)
print*,stat
stat=nf90_inq_varid(netid,'u',varid)
print*,stat
stat=nf90_get_var(netid,varid,ud,start=(/1,1,30,1/),count=(/nx,ny+1/))
print*,stat
stat=nf90_inq_varid(netid,'v',varid)
print*,stat
stat=nf90_get_var(netid,varid,vd,start=(/1,1,30,1/),count=(/nx+1,ny/))
print*,stat
stat=nf90_inq_varid(netid,'ua',varid)
print*,stat
stat=nf90_get_var(netid,varid,ua,start=(/1,1,30,1/),count=(/nx,ny/))
print*,stat
stat=nf90_inq_varid(netid,'va',varid)
print*,stat
stat=nf90_get_var(netid,varid,va,start=(/1,1,30,1/),count=(/nx,ny/))
print*,stat
print*,'lat=',maxval(lat2),minval(lat2)

! Pick off "B" grid points (corners of Grid cell)
DO j=1,nx+1
   DO i=1,nx+1
      lat(i,j)=lat2(2*i-1,2*j-1)
      lon(i,j)=lon2(2*i-1,2*j-1)
   ENDDO
ENDDO

! fudge u_w,v_w,u_s,v_s
u_s(:,1,:)=ua(:,1,:)
v_s(:,1,:)=va(:,1,:)
do j=2,ny
      u_s(:,j,:)=0.5*(ua(:,j-1,:)+ua(:,j,:))
      v_s(:,j,:)=0.5*(va(:,j-1,:)+va(:,j,:))
enddo
u_s(:,ny+1,:)=ua(:,ny,:)
v_s(:,ny+1,:)=va(:,ny,:)

u_w(1,:,:)=ua(1,:,:)
v_w(1,:,:)=va(1,:,:)
do i=2,nx
   u_w(i,:,:)=0.5*(ua(i-1,:,:)+ua(i,:,:))
   v_w(i,:,:)=0.5*(va(i-1,:,:)+va(i,:,:))
enddo
u_w(nx+1,:,:)=ua(nx,:,:)
v_w(nx+1,:,:)=va(nx,:,:)





! v point (left side)
do j=1,ny
   do i=1,nx+1
      p1(1) = lon(i,j)
      p1(2) = lat(i,j)
      p2(1) = lon(i,j+1)
      p2(2) = lat(i,j+1)

      call  mid_pt_sphere(p1, p2, p3)
      call get_unit_vect2(p1, p2, e1)
      call get_latlon_vector(p3, ex, ey)
      do k=1,km
         vdnew(i,j,k) = u_w(i,j,k)*inner_prod(e1,ex) + &
                        v_w(i,j,k)*inner_prod(e1,ey)
      enddo
   enddo
enddo



! upoint south
do j=1,ny+1
   do i=1,nx
      p1(1) = lon(i,j)
      p1(2) = lat(i,j)
      p2(1) = lon(i+1,j)
      p2(2) = lat(i+1,j)
      call  mid_pt_sphere(p1, p2, p3)
      call get_unit_vect2(p1, p2, e1)
      call get_latlon_vector(p3, ex, ey)
      do k=1,km
         udnew(i,j,k) = u_s(i,j,k)*inner_prod(e1,ex) + &
                        v_s(i,j,k)*inner_prod(e1,ey)
      enddo
   enddo
enddo
open(11,file='ud_tile'//tilef//'.bin',form='unformatted')
write(11) ud
write(11) udnew
write(11) u_s  
write(11) v_s  
close(11)
open(11,file='vd_tile'//tilef//'.bin',form='unformatted')
write(11) vd
write(11) vdnew
write(11) u_w  
write(11) v_w  
close(11)
open(11,file='awinds_tile'//tilef//'.bin',form='unformatted')
write(11) ua
write(11) va
close(11)


enddo
end

subroutine mid_pt_sphere(p1, p2, pm)
real , intent(IN)  :: p1(2), p2(2)
real , intent(OUT) :: pm(2)
real e1(3), e2(3), e3(3)

call latlon2xyz(p1, e1)
call latlon2xyz(p2, e2)
call mid_pt3_cart(e1, e2, e3)
call cart_to_latlon(1, e3, pm(1), pm(2))

end subroutine mid_pt_sphere

subroutine get_unit_vect2( e1, e2, uc )
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

real, intent(in) :: p(2)
real, intent(out):: e(3)


e(1) = cos(p(2)) * cos(p(1))
e(2) = cos(p(2)) * sin(p(1))
e(3) = sin(p(2))

end subroutine latlon2xyz
 
 
 
!>@brief The subroutine 'vect_cross performs cross products
!! of 3D vectors: e = P1 X P2
subroutine vect_cross(e, p1, p2)
real, intent(in) :: p1(3), p2(3)
real, intent(out):: e(3)

e(1) = p1(2)*p2(3) - p1(3)*p2(2)
e(2) = p1(3)*p2(1) - p1(1)*p2(3)
e(3) = p1(1)*p2(2) - p1(2)*p2(1)

end subroutine vect_cross


!>@brief The subroutine 'normalize_vect' makes 'e' a unit vector.
subroutine normalize_vect(e)

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
! vector version of cart_to_latlon1
integer, intent(in):: np
real, intent(inout):: q(3,np)
real, intent(inout):: xs(np), ys(np)
! local
real, parameter:: esl=1.d-10
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
real,intent(in):: v1(3), v2(3)
real :: prod16
integer k

prod16 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
inner_prod = prod16

end function inner_prod 
