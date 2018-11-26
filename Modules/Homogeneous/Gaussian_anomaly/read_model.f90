module read_model !!! Gaussian anomaly in a homogeneous medium !!! 


use angles

implicit none

public :: get_value, get_value_aniso!, get_value_anom, get_value_aniso_anom


contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,vp,vs,Qmu

doubleprecision :: r_c,theta_c,phi_c,theta_c_rad,phi_c_rad, lat,long, &
                   x_c,y_c,z_c, x_p,y_p,z_p, amp_anom,sigma,norm,anom


r_c = 6371000.d0 - 50000.d0   ! rayon du centre C de l'anomalie
theta_c = 91.d0   ! colatitude du centre C de l'anomalie
phi_c = -1.d0   ! longitude du centre C de l'anomalie
amp_anom = 0.1   ! amplitude de l'anomalie
sigma = 60000.d0   ! taille de l'anomalie en m (on aura zero a environ 3*sigma du centre)

! Coordonnees cartesiennes de C
theta_c_rad = deg2rad(theta_c)
phi_c_rad = deg2rad(phi_c)
x_c = r_c*dsin(theta_c_rad)*dcos(phi_c_rad)
y_c = r_c*dsin(theta_c_rad)*dsin(phi_c_rad)
z_c = r_c*dcos(theta_c_rad)

! Coordonnees cartesiennes du point d'entree P
x_p = r*dsin(theta_rad)*dcos(phi_rad)
y_p = r*dsin(theta_rad)*dsin(phi_rad)
z_p = r*dcos(theta_rad)

! Valeur de l'anomalie en P
norm = (x_p-x_c)**2 + (y_p-y_c)**2 + (z_p-z_c)**2
anom = amp_anom * dexp(-0.5*norm/sigma**2)

    vs = 5000.d0*(1.d0+anom)
    vp = 8000.d0
    rho = 3000.d0
    Qmu = 300.d0


end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu

doubleprecision :: vp,vs, vpv,vph,vsv,vsh, eta_aniso, &
                   r_c,theta_c,phi_c,theta_c_rad,phi_c_rad, lat,long, &
                   x_c,y_c,z_c, x_p,y_p,z_p, amp_anom,sigma,norm,anom


r_c = 6371000.d0 - 50000.d0   ! rayon du centre C de l'anomalie
theta_c = 91.d0   ! colatitude du centre C de l'anomalie
phi_c = -1.d0   ! longitude du centre C de l'anomalie
amp_anom = 0.1   ! amplitude de l'anomalie
sigma = 60000.d0   ! taille de l'anomalie en m (on aura zero a environ 3*sigma du centre)
    
! Coordonnees cartesiennes de C
theta_c_rad = deg2rad(theta_c)
phi_c_rad = deg2rad(phi_c)
x_c = r_c*dsin(theta_c_rad)*dcos(phi_c_rad)
y_c = r_c*dsin(theta_c_rad)*dsin(phi_c_rad)
z_c = r_c*dcos(theta_c_rad)
    
! Coordonnees cartesiennes du point d'entree P
x_p = r*dsin(theta_rad)*dcos(phi_rad)
y_p = r*dsin(theta_rad)*dsin(phi_rad)
z_p = r*dcos(theta_rad)

! Valeur de l'anomalie en P
norm = (x_p-x_c)**2 + (y_p-y_c)**2 + (z_p-z_c)**2
anom = amp_anom * dexp(-0.5*norm/sigma**2)

    vs = 5000.d0
    vp = 8000.d0
    rho = 3000.d0
    Qmu = 300.d0

    vsv = vs;   vsh = vs*(1.d0+anom)
    vpv = vp;   vph = vp
    eta_aniso = 1.d0

    A = rho*vph**2
    C = rho*vpv**2
    L = rho*vsv**2
    M = rho*vsh**2
    F = eta_aniso*(A-2.d0*L)
    Gc=0.d0; Gs=0.d0; Hc=0.d0; Hs=0.d0; Bc=0.d0; Bs=0.d0; Ec=0.d0; Es=0.d0


end subroutine get_value_aniso

!! ########################################################
!subroutine get_value_anom (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)
!
!
!implicit none
!
!integer, optional, intent(IN) :: moho
!doubleprecision, intent(IN) :: r,theta_rad,phi_rad
!doubleprecision, intent(OUT) :: rho,vp,vs,Qmu
!
!doubleprecision :: r_c,theta_c,phi_c,theta_c_rad,phi_c_rad, lat,long, &
!                   x_c,y_c,z_c, x_p,y_p,z_p, amp_anom,sigma,norm,anom
!
!
!r_c = 6371000.d0 - 50000.d0   ! rayon du centre C de l'anomalie
!theta_c = 91.d0   ! colatitude du centre C de l'anomalie
!phi_c = -1.d0   ! longitude du centre C de l'anomalie
!amp_anom = 0.1   ! amplitude de l'anomalie
!sigma = 60000.d0   ! taille de l'anomalie en m (on aura zero a environ 3*sigma du centre)
!
!! Coordonnees cartesiennes de C
!theta_c_rad = deg2rad(theta_c)
!phi_c_rad = deg2rad(phi_c)
!x_c = r_c*dsin(theta_c_rad)*dcos(phi_c_rad)
!y_c = r_c*dsin(theta_c_rad)*dsin(phi_c_rad)
!z_c = r_c*dcos(theta_c_rad)
!
!! Coordonnees cartesiennes du point d'entree P
!x_p = r*dsin(theta_rad)*dcos(phi_rad)
!y_p = r*dsin(theta_rad)*dsin(phi_rad)
!z_p = r*dcos(theta_rad)
!
!! Valeur de l'anomalie en P
!norm = (x_p-x_c)**2 + (y_p-y_c)**2 + (z_p-z_c)**2
!anom = amp_anom * dexp(-0.5*norm/sigma**2)
!
!    vs = 5000.d0*(1.d0+anom)
!    vp = 8000.d0
!    rho = 3000.d0
!    Qmu = 300.d0
!
!
!end subroutine get_value_anom
!
!! #######################################################
!subroutine get_value_aniso_anom (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)
!
!
!implicit none
!
!integer, optional, intent(IN) :: moho
!doubleprecision, intent(IN) :: r,theta_rad,phi_rad
!doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu
!
!doubleprecision :: vp,vs, vpv,vph,vsv,vsh, eta_aniso, &
!                   r_c,theta_c,phi_c,theta_c_rad,phi_c_rad, lat,long, &
!                   x_c,y_c,z_c, x_p,y_p,z_p, amp_anom,sigma,norm,anom
!
!
!r_c = 6371000.d0 - 50000.d0   ! rayon du centre C de l'anomalie
!theta_c = 91.d0   ! colatitude du centre C de l'anomalie
!phi_c = -1.d0   ! longitude du centre C de l'anomalie
!amp_anom = 0.1   ! amplitude de l'anomalie
!sigma = 60000.d0   ! taille de l'anomalie en m (on aura zero a environ 3*sigma du centre)
!
!! Coordonnees cartesiennes de C
!theta_c_rad = deg2rad(theta_c)
!phi_c_rad = deg2rad(phi_c)
!x_c = r_c*dsin(theta_c_rad)*dcos(phi_c_rad)
!y_c = r_c*dsin(theta_c_rad)*dsin(phi_c_rad)
!z_c = r_c*dcos(theta_c_rad)
!
!! Coordonnees cartesiennes du point d'entree P
!x_p = r*dsin(theta_rad)*dcos(phi_rad)
!y_p = r*dsin(theta_rad)*dsin(phi_rad)
!z_p = r*dcos(theta_rad)
!
!! Valeur de l'anomalie en P
!norm = (x_p-x_c)**2 + (y_p-y_c)**2 + (z_p-z_c)**2
!anom = amp_anom * dexp(-0.5*norm/sigma**2)
!
!
!    vs = 5000.d0
!    vp = 8000.d0
!    rho = 3000.d0
!    Qmu = 300.d0
!
!    vsv = vs;   vsh = vs*(1.d0+anom)
!    vpv = vp;   vph = vp
!    eta_aniso = 1.d0
!
!    A = rho*vph**2
!    C = rho*vpv**2
!    L = rho*vsv**2
!    M = rho*vsh**2
!    F = eta_aniso*(A-2.d0*L)
!    Gc=0.d0; Gs=0.d0; Hc=0.d0; Hs=0.d0; Bc=0.d0; Bs=0.d0; Ec=0.d0; Es=0.d0
!
!
!end subroutine get_value_aniso_anom

! #######################################################
end module read_model
