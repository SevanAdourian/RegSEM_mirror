module read_model !!! Two homogeneous layers !!! 


use angles

implicit none

public :: get_value, get_value_aniso


contains

! #######################################################
subroutine get_value (r,theta_rad,phi_rad,rho,vp,vs,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,vp,vs,Qmu

doubleprecision, parameter :: interface_radius = 6366000.d0

if (r>interface_radius) then
    vs = 3850.d0
    vp = 6750.d0
    rho = 3000.d0
    Qmu = 300.d0
else
    vs = 4400.d0
    vp = 8000.d0
    rho = 3400.d0
    Qmu = 200.d0
endif


end subroutine get_value

! #######################################################
subroutine get_value_aniso (r,theta_rad,phi_rad,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu,moho)


implicit none

integer, optional, intent(IN) :: moho
doubleprecision, intent(IN) :: r,theta_rad,phi_rad
doubleprecision, intent(OUT) :: rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu

doubleprecision :: vp,vs, vsv,vsh,vpv,vph, eta_aniso, A0,C0,F0,L0,M0
doubleprecision, parameter :: interface_radius = 6341000.d0, &
!doubleprecision, parameter :: interface_radius = 17000.d0, &
                              anom_vs = 0.2, anom_vp = 0.2

if (r>interface_radius) then
!if (phi_rad>interface_radius) then   ! Dans un cube, phi_rad est la coordonnee z (i.e. la verticale)
    vs = 3850.d0;   vsv = vs*(1.d0+anom_vs);   vsh = vs*(1.d0-anom_vs)
    vp = 6750.d0;   vpv = vp*(1.d0+anom_vp);   vph = vp*(1.d0-anom_vp)
    eta_aniso = 1.1
    rho = 3000.d0
    Qmu = 300.d0
else
    vs = 4400.d0;   vsv = vs;   vsh = vs 
    vp = 8000.d0;   vpv = vp;   vph = vp
    eta_aniso = 1.d0
    rho = 3400.d0
    Qmu = 200.d0
!    vs = 3850.d0;   vsv = vs;   vsh = vs
!    vp = 6750.d0;   vpv = vp;   vph = vp
!    eta_aniso = 1.d0
!    rho = 3000.d0
!    Qmu = 300.d0
endif 

!A = rho*vph**2
!C = rho*vpv**2
!L = rho*vsv**2
!M = rho*vsh**2
!F = eta_aniso*(A-2.d0*L)
!Gc=0.d0; Gs=0.d0; Hc=0.d0; Hs=0.d0; Bc=0.d0; Bs=0.d0; Ec=0.d0; Es=0.d0

A0 = rho*vph**2
C0 = rho*vpv**2
L0 = rho*vsv**2
M0 = rho*vsh**2
F0 = eta_aniso*(A0-2.d0*L0)

! VTI -> HTI est-ouest
A = 3.d0*(A0+C0)/8.d0 + F0/4.d0 + L0/2.d0
C = A0
F = (A0+F0)/2.d0 - M0
L = (M0+L0)/2.d0
M = (A0+C0)/8.d0 - F0/4.d0 + L0/2.d0
Gc = (M0-L0)/2.d0; Gs = 0.d0
Hc = (A0-F0)/2.d0 - M0; Hs = 0.d0
Bc = (A0-C0)/2.d0; Bs = 0.d0
Ec = (A0+C0)/8.d0 - F0/4.d0 - L0/2.d0; Es = 0.d0


end subroutine get_value_aniso

! #######################################################
end module read_model
