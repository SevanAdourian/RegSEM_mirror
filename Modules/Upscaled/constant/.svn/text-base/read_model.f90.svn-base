module read_model !!! For a model obtained from homogenization of a periodic medium !!!
! Routine get_homog_param is used here instead of get_value_aniso.
! In define_arrays.f90, lines 279 to 298 have to be commented.

implicit none

public :: get_value, get_value_aniso, get_homog_param
private :: rho_homog, C_homog, first_time

doubleprecision :: rho_homog
doubleprecision, dimension(1:21) :: C_homog
logical :: first_time = .true.

contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,vp,vs,Qmu


stop 'THE CUBES ARE ANISOTROPIC !!!'


end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)

implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu


stop 'PLEASE USE get_homog_param INSTEAD OF get_value_aniso !!!'


end subroutine get_value_aniso

! ########################################################
subroutine get_homog_param (rho,C)

implicit none

doubleprecision, intent(OUT) :: rho
doubleprecision, dimension(1:6,1:6), intent(OUT) :: C

integer :: i,j, ios
doubleprecision, dimension(1:2) :: tmp


if (first_time) then
   first_time = .false.
   open (11,file="homog_param.txt",form="formatted",status="old",iostat=ios)
   if (ios>0)   stop 'TEXT FILE homog_param.txt IS NEEDED TO GET THE ELASTIC PARAMETERS OF THE MEDIUM'
   do i = 1,21
      read (11,*) C_homog(i), tmp(1:2)   ! NB: Cij values in homog_param.txt are in Mandel's notation
   enddo
   read (11,*) rho_homog, tmp(1:2)
   close (11)
endif

rho = rho_homog
C(1,1:6) = C_homog(1:6)
C(2,2:6) = C_homog(7:11)
C(3,3:6) = C_homog(12:15)
C(4,4:6) = C_homog(16:18);
C(5,5:6) = C_homog(19:20);
C(6,6) = C_homog(21);
do i = 2,6
   do j = 1,i-1
      C(i,j) = C(j,i)
   enddo
enddo


end subroutine get_homog_param

! #######################################################
end module read_model
