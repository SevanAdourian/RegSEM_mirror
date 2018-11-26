module read_model !!! Polynomial interpolation of a homogenized non-periodic medium!!!
! Routine get_elastic_param is used here instead of get_value_aniso.
! In define_arrays.f90, lines 279 to 298 have to be commented.

implicit none

public :: get_value, get_value_aniso, get_elastic_param, deload
private :: nb_pt_dir, model_bounds, rho_homog, C_homog, xseries,yseries,zseries, &
           first_time, order

integer, parameter :: order = 4   ! 2=linear, 3=quadratic, 4=cubic, etc...

integer, dimension(1:3) :: nb_pt_dir
doubleprecision, dimension(:), allocatable :: xseries,yseries,zseries
doubleprecision, dimension(1:3,1:2) :: model_bounds
doubleprecision, dimension(:,:,:), allocatable :: rho_homog
doubleprecision, dimension(:,:,:,:), allocatable :: C_homog
logical :: first_time = .true.

contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,vp,vs,Qmu


stop 'THE UPSCALED MEDIUM IS ANISOTROPIC !!!'


end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu


stop 'PLEASE USE get_elastic_param INSTEAD OF get_value_aniso !!!'


end subroutine get_value_aniso

! ########################################################
subroutine get_elastic_param (x,y,z,rho,C,dom_bounds)

use module_polinterp

implicit none

doubleprecision, intent(IN) :: x,y,z
doubleprecision, intent(OUT) :: rho
doubleprecision, dimension(1:3,1:2), intent(IN) :: dom_bounds
doubleprecision, dimension(1:6,1:6), intent(OUT) :: C

integer :: i,j,k, ii,jj,kk, ios, i_count, length
integer, dimension(1:3,1:2) :: limit
doubleprecision :: debut, lambda, mu
doubleprecision, dimension(:), allocatable :: tmp
doubleprecision, dimension(1:3) :: delta_dir
doubleprecision, dimension(1:21) :: tmp_C

integer, parameter :: nb_pt_val = 22   ! rho + C


if (first_time) then

   first_time = .false.
   allocate (tmp(1:nb_pt_val))
   inquire(iolength=length) tmp
   open (11,file="Cstar_rho.res",form="unformatted",access="direct",recl=length,status="old",iostat=ios)
   if (ios>0)   stop 'BINARY FILE Cstar_rho.res IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'

   read (11,rec=1) tmp
   i_count = 1
   do i = 1,3
      model_bounds(i,1) = tmp(i_count);   model_bounds(i,2) = tmp(i_count+1);   nb_pt_dir(i) = int(tmp(i_count+2))
      delta_dir(i) = (model_bounds(i,2)-model_bounds(i,1)) / (nb_pt_dir(i)-1)
      limit(i,1) = int( (dom_bounds(i,1)-model_bounds(i,1)) / delta_dir(i) ) + 1
      if (limit(i,1)<1)   limit(i,1) = 1
      if (limit(i,1)>nb_pt_dir(i))   limit(i,1) = nb_pt_dir(i)
      limit(i,2) = int( (dom_bounds(i,2)-model_bounds(i,1)) / delta_dir(i) ) + 2
      if (limit(i,2)<1)   limit(i,2) = 1
      if (limit(i,2)>nb_pt_dir(i))   limit(i,2) = nb_pt_dir(i)
      i_count = i_count + 3
   enddo

   ii = limit(1,2)-limit(1,1)+1
   jj = limit(2,2)-limit(2,1)+1
   kk = limit(3,2)-limit(3,1)+1
   allocate (rho_homog(1:ii,1:jj,1:kk))
   allocate (C_homog(1:ii,1:jj,1:kk,1:21))
   do k = limit(3,1),limit(3,2)
      kk = k - limit(3,1) + 1
      do j = limit(2,1),limit(2,2)
         jj = j - limit(2,1) + 1
         do i = limit(1,1),limit(1,2)
            ii = i - limit(1,1) + 1
            i_count = (k-1)*nb_pt_dir(1)*nb_pt_dir(2) + &
                      (j-1)*nb_pt_dir(1) + i + 1   ! +1 because of the header
            read (11,rec=i_count) tmp
            rho_homog(ii,jj,kk) = tmp(1)
            C_homog(ii,jj,kk,1:21) = tmp(2:22)
!lambda = tmp(2) + tmp(8) + tmp(13) + 4.d0*(tmp(3)+tmp(4)+tmp(9)) - (tmp(17)+tmp(20)+tmp(22))
!lambda = lambda/15.d0
!mu = tmp(2) + tmp(8) + tmp(13) - (tmp(3)+tmp(4)+tmp(9)) + 1.5*(tmp(17)+tmp(20)+tmp(22))
!mu = mu/15.d0
!if (dsqrt((lambda+2*mu)/tmp(1))>20000.d0) then
!   print *, "SHIT!!! : ", dsqrt((lambda+2*mu)/tmp(1))
!   print *, i,j,k
!endif
         enddo
      enddo
   enddo
   close (11)

   ii = limit(1,2)-limit(1,1)+1
   jj = limit(2,2)-limit(2,1)+1
   kk = limit(3,2)-limit(3,1)+1
   allocate (xseries(1:ii))
   allocate (yseries(1:jj))
   allocate (zseries(1:kk))
   debut = model_bounds(1,1) + (limit(1,1)-1)*delta_dir(1)
   do i = 1,ii
      xseries(i) = debut + (i-1)*delta_dir(1)
   enddo
   debut = model_bounds(2,1) + (limit(2,1)-1)*delta_dir(2)
   do j = 1,jj
      yseries(j) = debut + (j-1)*delta_dir(2)
   enddo
   debut = model_bounds(3,1) + (limit(3,1)-1)*delta_dir(3)
   do k = 1,kk
      zseries(k) = debut + (k-1)*delta_dir(3)
   enddo

endif

if (x<model_bounds(1,1) .or. x>model_bounds(1,2) .or. &
    y<model_bounds(2,1) .or. y>model_bounds(2,2) .or. &
    z<model_bounds(3,1) .or. z>model_bounds(3,2)) then
   rho = 3003.74107037824;   lambda = 55260199143.1541;   mu = 42994835872.9560
   tmp_C(:) = 0.d0
   tmp_C(1) = 2.d0*mu + lambda;   tmp_C(7) = 2.d0*mu + lambda;   tmp_C(12) = 2.d0*mu + lambda
   tmp_C(16) = 2.d0*mu;   tmp_C(19) = 2.d0*mu;   tmp_C(21) = 2.d0*mu
   tmp_C(2) = lambda;   tmp_C(3) = lambda;   tmp_C(8) = lambda
else
   do i_count = 1,21
      call polint3D (xseries, yseries, zseries, C_homog(:,:,:,i_count), order, x, y, z, tmp_C(i_count))
   enddo
   call polint3D (xseries, yseries, zseries, rho_homog(:,:,:), order, x, y, z, rho)
endif

C(1,1:6) = tmp_C(1:6)
C(2,2:6) = tmp_C(7:11)
C(3,3:6) = tmp_C(12:15)
C(4,4:6) = tmp_C(16:18)
C(5,5:6) = tmp_C(19:20)
C(6,6) = tmp_C(21)
do i = 2,6
   do j = 1,i-1
      C(i,j) = C(j,i)
   enddo
enddo
!lambda = tmp_C(1) + tmp_C(7) + tmp_C(12) + 4.d0*(tmp_C(2)+tmp_C(3)+tmp_C(8)) - (tmp_C(16)+tmp_C(19)+tmp_C(21))
!lambda = lambda/15.d0
!mu = tmp_C(1) + tmp_C(7) + tmp_C(12) - (tmp_C(2)+tmp_C(3)+tmp_C(8)) + 1.5*(tmp_C(16)+tmp_C(19)+tmp_C(21))
!mu = mu/15.d0
!if (dsqrt((lambda+2*mu)/rho)>20000.d0) then
!   print *, "AIE!!! : ", dsqrt((lambda+2*mu)/rho)
!   print *, x,y,z
!endif


end subroutine get_elastic_param

! #######################################################
subroutine deload

implicit none


deallocate (rho_homog, C_homog, xseries, yseries, zseries)


end subroutine deload


end module read_model
