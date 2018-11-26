module read_model !!! Spline interpolation of a homogenized non-periodic medium!!!
! Routine get_elastic_param is used here instead of get_value_aniso.
! In define_arrays.f90, lines 279 to 298 have to be commented.

implicit none

public :: get_value, get_value_aniso, get_elastic_param
private :: nb_pt_dir, rho_homog,rho_spline, C_homog,C_spline, xseries,yseries,zseries, &
           y_rho_val,z_rho_val, y_spline,z_spline, y_C_val,z_C_val, first_time

integer, dimension(1:3) :: nb_pt_dir
doubleprecision, dimension(:), allocatable :: xseries,yseries,zseries, &
                                              y_rho_val,z_rho_val, y_spline,z_spline
doubleprecision, dimension(:,:), allocatable :: y_C_val,z_C_val
doubleprecision, dimension(:,:,:), allocatable :: rho_homog, rho_spline
doubleprecision, dimension(:,:,:,:), allocatable :: C_homog, C_spline
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
subroutine get_elastic_param (x,y,z,rho,dom_bounds,C)

use module_spline

implicit none

doubleprecision, intent(IN) :: x,y,z
doubleprecision, intent(OUT) :: rho
doubleprecision, dimension(1:3,1:2), intent(IN) :: dom_bounds   ! not used yet
doubleprecision, dimension(1:6,1:6), intent(OUT) :: C

integer :: i,j,k, ios, i_count
doubleprecision, dimension(1:3) :: delta_dir, tmp
doubleprecision, dimension(1:21) :: tmp_C
doubleprecision, dimension(1:3,1:2) :: bounds


if (first_time) then
   first_time = .false.
   !!! The following won't work with the new format of Cstar_rho.res (direct access)
   open (11,file="Cstar_rho.res",form="unformatted",access="sequential",status="old",iostat=ios)
   if (ios>0)   stop 'TEXT FILE Cstar_rho.res IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'
   do i = 1,3
      read (11) bounds(i,1), bounds(i,2), nb_pt_dir(i)
      delta_dir(i) = (bounds(i,2)-bounds(i,1)) / (nb_pt_dir(i)-1)
   enddo
   allocate (C_homog(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3),1:21))
   allocate (rho_homog(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3)))
   do k = 1,nb_pt_dir(3)
    do j = 1,nb_pt_dir(2)
     do i = 1,nb_pt_dir(1)
        read (11) tmp(1:3)
        read (11) rho_homog(i,j,k)
        do i_count = 1,21
           read (11) C_homog(i,j,k,i_count)   ! NB: Cij values in Cstar_rho.res are in Mandel's notation
        enddo
        read (11) tmp(1)
     enddo
    enddo
   enddo
   close (11)
   allocate (xseries(1:nb_pt_dir(1)))
   allocate (yseries(1:nb_pt_dir(2)))
   allocate (zseries(1:nb_pt_dir(3)))
   do i = 1,nb_pt_dir(1)
      xseries(i) = bounds(1,1) + (i-1)*delta_dir(1)
   enddo
   do i = 1,nb_pt_dir(2)
      yseries(i) = bounds(2,1) + (i-1)*delta_dir(2)
   enddo
   do i = 1,nb_pt_dir(3)
      zseries(i) = bounds(3,1) + (i-1)*delta_dir(3)
   enddo
   allocate (C_spline(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3),1:21))
   allocate (rho_spline(1:nb_pt_dir(1),1:nb_pt_dir(2),1:nb_pt_dir(3)))
   do k = 1,nb_pt_dir(3)
      do j = 1,nb_pt_dir(2)
         do i_count = 1,21
            call spline (xseries,C_homog(:,j,k,i_count),1.d30,1.d30,C_spline(:,j,k,i_count))
         enddo
         call spline (xseries,rho_homog(:,j,k),1.d30,1.d30,rho_spline(:,j,k))
      enddo
   enddo
   allocate (y_C_val(1:nb_pt_dir(2),1:21))
   allocate (y_rho_val(1:nb_pt_dir(2)))
   allocate (y_spline(1:nb_pt_dir(2)))
   allocate (z_C_val(1:nb_pt_dir(3),1:21))
   allocate (z_rho_val(1:nb_pt_dir(3)))
   allocate (z_spline(1:nb_pt_dir(3)))
endif

do k = 1,nb_pt_dir(3)
   do j = 1,nb_pt_dir(2)
      do i_count = 1,21
         y_C_val(j,i_count) = splintx(xseries,C_homog(:,j,k,i_count),C_spline(:,j,k,i_count),x)
      enddo
      y_rho_val(j) = splintx(xseries,rho_homog(:,j,k),rho_spline(:,j,k),x)
   enddo
   do i_count = 1,21
      call spline (yseries,y_C_val(:,i_count),1.d30,1.d30,y_spline)
      z_C_val(k,i_count) = splinty(yseries,y_C_val(:,i_count),y_spline,y)
   enddo
   call spline (yseries,y_rho_val,1.d30,1.d30,y_spline)
   z_rho_val(k) = splinty(yseries,y_rho_val,y_spline,y)
enddo
do i_count = 1,21
   call spline (zseries,z_C_val(:,i_count),1.d30,1.d30,z_spline)
   tmp_C(i_count) = splintz(zseries,z_C_val(:,i_count),z_spline,z)
enddo
call spline (zseries,z_rho_val,1.d30,1.d30,z_spline)
rho = splintz(zseries,z_rho_val,z_spline,z)

C(1,1:6) = tmp_C(1:6)
C(2,2:6) = tmp_C(7:11)
C(3,3:6) = tmp_C(12:15)
C(4,4:6) = tmp_C(16:18);
C(5,5:6) = tmp_C(19:20);
C(6,6) = tmp_C(21);
do i = 2,6
   do j = 1,i-1
      C(i,j) = C(j,i)
   enddo
enddo


end subroutine get_elastic_param

! #######################################################
end module read_model
