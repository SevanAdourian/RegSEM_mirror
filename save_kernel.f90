subroutine save_kernel (Tdomain,rg,delta_chi)


use sdomains
use angles

implicit none

include 'mpif.h'

type (Domain), intent (INOUT) :: Tdomain
integer, intent(IN) :: rg
doubleprecision, intent(INOUT) :: delta_chi

!! ------------------------------------------------------------------
!! compile-time constant: include scaled perturbations for TI medium?
logical, parameter :: use_ma89_scalings = .true.
!! ------------------------------------------------------------------

integer :: n, x,y,z, ngll1,ngll2,ngll3, i,j,k, code
integer :: mat !! sfrench
doubleprecision :: u,v,w, xa,ya,za, r,theta,phi, lat_deg,phi_deg, delta, &
                   rho,lambda,mu, A,C,L,M, kern_beta,kern_xi, &
                   kern_alpha,kern_phi, kern_eta, kern_rho, KB, KX
doubleprecision :: wgt !! sfrench
doubleprecision, dimension(0:2) :: tmp
doubleprecision, dimension (0:2,0:2) :: Rot
doubleprecision, allocatable, dimension(:,:,:) :: Rsph !! sfrench
character*60 :: fname1, fname2, name_vtk_file
logical, parameter :: write_vtk_kernels = .true.
integer :: imin, imax, jmin, jmax, kmin, kmax, indx, indy, indz, nrec
INTEGER X1, X2, Y1, Y2, Z1, Z2, I1, I2, J1, J2, K1, K2, IPROC
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XC,YC,ZC
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: kernel_VALUES
integer :: E_IO
integer, parameter :: Unit_VTK=27
integer ::  offset
integer :: dim_arrays
character(1), parameter:: end_rec = char(10)
character*20 :: cx1, cx2, cy1, cy2, cz1, cz2, coffset
double precision :: scale_beta,scale_xi
integer :: indxe, indye, indze, i1e, i2e, j1e, j2e, k1e, k2e

i1e = Tdomain%i1_mirror
i2e = Tdomain%n_elem_x-Tdomain%i2_mirror-1
j1e = Tdomain%j1_mirror
j2e = Tdomain%n_elem_y-Tdomain%j2_mirror-1
k1e = Tdomain%k1_mirror
k2e = Tdomain%n_elem_z-Tdomain%k2_mirror-1

write (fname1,"(a,I3.3)") "kernel_beta_",rg
open (31,file=trim(fname1))
if (Tdomain%aniso) then
    write (fname2,"(a,I3.3)") "kernel_xi_",rg
    open (32,file=trim(fname2))
endif

! for paraview output

if(write_vtk_kernels)then
     write (name_vtk_file,"(a,I3.3)") "kernels_for_vtk_",rg
     open (26,file=trim(name_vtk_file),form='UNFORMATTED',access='stream')
     do i = 1,7
        write(26)nrec
     enddo
     
     nrec=0
     
     imin =  999999.
     imax = -999999.
     jmin =  999999.
     jmax = -999999.
     kmin =  999999.
     kmax = -999999.

endif

delta_chi = 0.d0
scale_beta = 0.d0
scale_xi = 0.d0

do n = 0,Tdomain%n_elem-1

    indxe = Tdomain%specel(n)%indx
    indye = Tdomain%specel(n)%indy
    indze = Tdomain%specel(n)%indz

    if (.not. Tdomain%specel(n)%PML) then

        ngll1 = Tdomain%specel(n)%ngllx
        ngll2 = Tdomain%specel(n)%nglly
        ngll3 = Tdomain%specel(n)%ngllz

        !! sfrench
        if (Tdomain%ellipticity) then
            ! fetch spherical-earth radius as well if have ellipticity
            allocate (Rsph(0:ngll1-1,0:ngll2-1,0:ngll3-1))
            mat = Tdomain%specel(n)%mat_index
            call r_spheric(ngll1,ngll2,ngll3,Tdomain,mat,n,Rsph)
        endif

        indz = Tdomain%indz_gll(Tdomain%specel(n)%indz)+Tdomain%specel(n)%indz
        do z = 0,ngll3-1 
         indy = Tdomain%indy_gll(Tdomain%specel(n)%indy)+Tdomain%specel(n)%indy  
         do y = 0,ngll2-1 
          indx = Tdomain%indx_gll(Tdomain%specel(n)%indx)+Tdomain%specel(n)%indx
          do x = 0,ngll1-1
  
              nrec = nrec+1

              kmin = min(kmin,indz)
              kmax = max(kmax,indz)
              jmin = min(jmin,indy)
              jmax = max(jmax,indy)
              imin = min(imin,indx)
              imax = max(imax,indx)

              ! Taking the cartesian coordinates of the GLL
              i = Tdomain%specel(n)%Iglobnum(x,y,z)
              u = Tdomain%Globcoord(0,i)
              v = Tdomain%Globcoord(1,i)
              w = Tdomain%Globcoord(2,i)
              if (Tdomain%curve) then
                  ! Coordinates in the real chunk
                  xa = u;   ya = v;   za = w
                  Rot = Tdomain%rot
                  u = Rot(0,0)*xa + Rot(0,1)*ya + Rot(0,2)*za
                  v = Rot(1,0)*xa + Rot(1,1)*ya + Rot(1,2)*za
                  w = Rot(2,0)*xa + Rot(2,1)*ya + Rot(2,2)*za
                  ! Convert the cartesian coordinates into spherical coordinates
                  call cart2sph (u,v,w,r,theta,phi)
                  lat_deg = 90.d0 - 180.d0*theta/pi
                  phi_deg = 180.d0*phi/pi
                  if (phi_deg>180.d0)   phi_deg = phi_deg - 360.d0
                  tmp(0) = lat_deg;   tmp(1) = phi_deg;   tmp(2) = r/1000.d0
              else
                  tmp(0) = u/1000.d0;   tmp(1) = v/1000.d0;   tmp(2) = w/1000.d0
              endif

              
              if(write_vtk_kernels)write(26) indx, indy, indz, sngl(u), sngl(v), sngl(w)        

              !!! Computing beta and xi kernels !!!
              if (Tdomain%aniso) then
                  A = Tdomain%specel(n)%save_TIparam(1,x,y,z)
                  C = Tdomain%specel(n)%save_TIparam(2,x,y,z)
                  L = Tdomain%specel(n)%save_TIparam(4,x,y,z)
                  M = Tdomain%specel(n)%save_TIparam(5,x,y,z)
                  ! Matches Panning and Romanowicz (2006) Appendix A
                  ! eq. A.16
                  !write(*,*) "kern4 is", Tdomain%specel(n)%Kern_aniso(4,x,y,z)
                  !write(*,*) "kern5 is", Tdomain%specel(n)%Kern_aniso(5,x,y,z)
                  !write(*,*) "kern3 is", Tdomain%specel(n)%Kern_aniso(3,x,y,z)
                  !write(*,*) "A is", A
                  !write(*,*) "L is", L
                  !write(*,*) "C is", C
                  !write(*,*) "M is", M
                  kern_beta = 2.d0*( &
                            Tdomain%specel(n)%Kern_aniso(4,x,y,z) + &
                            Tdomain%specel(n)%Kern_aniso(5,x,y,z) - &
                      2.d0*L*Tdomain%specel(n)%Kern_aniso(3,x,y,z)/(A-2.d0*L))
                  ! eq. A.18
                  kern_xi   = ( &
                        2.d0*L*Tdomain%specel(n)%Kern_aniso(5,x,y,z) - &
                             M*Tdomain%specel(n)%Kern_aniso(4,x,y,z) + &
                      2.d0*L*M*Tdomain%specel(n)%Kern_aniso(3,x,y,z)/(A-2.d0*L))/(2.d0*L+M)
                  !
                  ! Possibly apply scalings (Montagner and Anderson, 1989)
                  if (use_ma89_scalings) then
                      ! eq. A.17
                      kern_alpha = 2.d0*( &
                            Tdomain%specel(n)%Kern_aniso(1,x,y,z) + &
                            Tdomain%specel(n)%Kern_aniso(2,x,y,z) + &
                          A*Tdomain%specel(n)%Kern_aniso(3,x,y,z)/(A-2.d0*L))
                      ! eq. A.19
                      kern_phi   = ( &
                          4.d0*A*Tdomain%specel(n)%Kern_aniso(2,x,y,z) - &
                               C*Tdomain%specel(n)%Kern_aniso(1,x,y,z) - &
                             A*C*Tdomain%specel(n)%Kern_aniso(3,x,y,z)/(A-2.d0*L))/(4.d0*A+C)
                      ! eq. A.20
                      kern_eta   = Tdomain%specel(n)%Kern_aniso(3,x,y,z)
                      ! eq. A.21
                      kern_rho   = &
                          Tdomain%specel(n)%Kern_rho(x,y,z) + &
                          Tdomain%specel(n)%Kern_aniso(1,x,y,z) + &
                          Tdomain%specel(n)%Kern_aniso(2,x,y,z) + &
                          Tdomain%specel(n)%Kern_aniso(3,x,y,z) + &
                          Tdomain%specel(n)%Kern_aniso(4,x,y,z) + &
                          Tdomain%specel(n)%Kern_aniso(5,x,y,z)
                      !
                      ! now add in the scaled kernels
                      kern_beta = kern_beta   + &
                          0.50d0 * kern_alpha + &
                          0.33d0 * kern_rho
                      kern_xi = kern_xi     - &
                          2.50d0 * kern_eta - &
                          1.50d0 * kern_phi
                  endif
              else
                  lambda = Tdomain%specel(n)%Lambda(x,y,z)
                  mu = Tdomain%specel(n)%Mu(x,y,z)
                  kern_beta = 2.d0 * (Tdomain%specel(n)%Kern_mu(x,y,z) - &
                                      2.d0*mu*Tdomain%specel(n)%Kern_lambda(x,y,z)/lambda)
              endif

              !!! check de la derivee partielle !!!
!              if (Tdomain%aniso) then
!                  delta = Tdomain%specel(n)%anomaly(x,y,z)/(M/L) - 1.d0
!                  delta_chi = delta_chi + &
!                              delta * kern_xi * &
!                              Tdomain%specel(n)%Jacob(x,y,z) * &
!                              Tdomain%specel(n)%wgtx(x) * Tdomain%specel(n)%wgty(y) * Tdomain%specel(n)%wgtz(z)
!              else
!                  rho = Tdomain%specel(n)%Density(x,y,z) 
!                  delta = Tdomain%specel(n)%anomaly(x,y,z)/dsqrt(mu/rho) - 1.d0
!                  delta_chi = delta_chi + &
!                              delta * kern_beta * &
!                              Tdomain%specel(n)%Jacob(x,y,z) * &
!                              Tdomain%specel(n)%wgtx(x) * Tdomain%specel(n)%wgty(y) * Tdomain%specel(n)%wgtz(z)
!              endif

              !!!!!! Storing the value of KERNELthe kernel at the interior GLLs !!!
              !!!if (x/=0 .and. y/=0 .and. z/=0 .and. x/=ngll1-1 .and. y/=ngll2-1 .and. z/=ngll3-1) then
              !!!    write (31,*) tmp(0:2),Tdomain%specel(n)%Jacob(x,y,z), kern_beta
              !!!    if (Tdomain%aniso)   write (32,*) tmp(0:2),Tdomain%specel(n)%Jacob(x,y,z), kern_xi
              !!!endif

              !! sfrench: now storing the value at _all_ GLLs for easy
              !!          projection onto model basis
              wgt = Tdomain%specel(n)%Jacob(x,y,z) * &
                    Tdomain%specel(n)%wgtx(x) * &
                    Tdomain%specel(n)%wgty(y) * &
                    Tdomain%specel(n)%wgtz(z)
              !! zero kernels outside window domain

              if(Tdomain%specel(n)%mirror_position<=1)then
               kern_beta = 0.d0
               kern_xi = 0.d0
              endif

              ! quick fix for grqdient borders
                ! ADD BY PIERRE BORDERS ON X AND Y are still present
                ! Was set to 1 and on Z to 0
              if(indxe<=i1e+2)kern_beta = 0.d0; kern_xi = 0.d0
              if(indxe>=i2e-2)kern_beta = 0.d0; kern_xi = 0.d0
              if(indye<=j1e+2)kern_beta = 0.d0; kern_xi = 0.d0
              if(indye>=j2e-2)kern_beta = 0.d0; kern_xi = 0.d0
              if(indze<=k1e+0)kern_beta = 0.d0; kern_xi = 0.d0
              if(indze>=k2e-0)kern_beta = 0.d0; kern_xi = 0.d0

              !! fmt: elem,n1,n2,n3,lat,lon,r,[r_sphere],weights * Jac,knl
              if (Tdomain%ellipticity) then
                  write (31,'(i5,x,3(i2,x),6(e14.4,x))') n,x,y,z,tmp(0:2),Rsph(x,y,z),wgt,kern_beta
                  if (Tdomain%aniso) &
                      write (32,'(i5,x,3(i2,x),6(e14.4,x))') n,x,y,z,tmp(0:2),Rsph(x,y,z),wgt,kern_xi
              else
                  write (31,'(i5,x,3(i2,x),5(e14.4,x))') n,x,y,z,tmp(0:2),wgt,kern_beta
                  if (Tdomain%aniso) &
                      write (32,'(i5,x,3(i2,x),5(e14.4,x))') n,x,y,z,tmp(0:2),wgt,kern_xi
              endif
          scale_beta = max(scale_beta,abs(kern_beta))
          scale_xi = max(scale_xi,abs(kern_xi))
          if(write_vtk_kernels)write(26)kern_beta, kern_xi

          indx = indx+1
          enddo
         indy = indy+1
         enddo
        indz = indz+1
        enddo

        !! sfrench
        if (Tdomain%ellipticity) deallocate(Rsph)

    endif

enddo

if(write_vtk_kernels)then 
  write(26,pos=1)nrec,imin,imax,jmin,jmax,kmin,kmax
  close (26)
  call MPI_BARRIER(MPI_COMM_WORLD, code)       
  if(rg==0)then
  !
  ! READ ALL FILE HEADER AND DETERMINE GRID DIMENSIONS
  !
  DO IPROC = 0,Tdomain%n_proc-1
    write (name_vtk_file,"(a,I3.3)") "kernels_for_vtk_",iproc
    open (26,file=trim(name_vtk_file),status='old',form='UNFORMATTED',access='stream')
     READ(26)NREC,I1,I2,J1,J2,K1,K2
     IF(IPROC==0)THEN
        X1=I1 ; X2=I2 ; Y1=J1 ; Y2=J2 ; Z1=K1 ; Z2=K2 
     ELSE
        X1 = MIN(X1,I1)
        X2 = MAX(X2,I2)
        Y1 = MIN(Y1,J1)
        Y2 = MAX(Y2,J2)
        Z1 = MIN(Z1,K1)
        Z2 = MAX(Z2,K2)
     ENDIF
     CLOSE(26)
  ENDDO
  !
  ALLOCATE(XC(x1:x2,y1:y2,z1:z2))
  ALLOCATE(YC(x1:x2,y1:y2,z1:z2))
  ALLOCATE(ZC(x1:x2,y1:y2,z1:z2))
  !
  ALLOCATE(KERNEL_VALUES(0:1,x1:x2,y1:y2,z1:z2))
  !
  ! MERGE FILES, READ/WRITE DATA
  !
  DO IPROC = 0,Tdomain%n_proc-1 
     write (NAME_VTK_FILE,"(a,I3.3)") "kernels_for_vtk_",iproc
     open (26,file=trim(NAME_VTK_FILE),status='old',form='UNFORMATTED',access='stream')
     READ(26)NREC,I1,I2,J1,J2,K1,K2     
     !      
     DO N = 1,NREC
        !
        READ(26)I        , J        , K
        READ(26)XC(I,J,K) , YC(I,J,K) , ZC(I,J,K)
        !READ(26)KERNEL_VALUES(:,i,j,k)
        READ(26)KB,KX
        KERNEL_VALUES(0,I,J,K) = SNGL(KB)
        KERNEL_VALUES(1,I,J,K) = SNGL(KX)
        !         
     ENDDO
     !      
     CLOSE(26, STATUS='DELETE')
     !
  ENDDO
  !
  ! METERS TO KILOMETERS
  ! 
  xC = XC/6371000.d0
  yC = YC/6371000.d0
  zC = ZC/6371000.d0
!
write(cx1,'(I10)')x1
write(cx2,'(I10)')x2
write(cy1,'(I10)')y1
write(cy2,'(I10)')y2
write(cz1,'(I10)')z1
write(cz2,'(I10)')z2
!
open(unit       = Unit_VTK,        &
     file       = 'kernels.vts',   &
     form       = 'UNFORMATTED',   &
     access     = 'STREAM',        &
     action     = 'WRITE',         &
     status     = 'REPLACE',       &
     convert    = 'LITTLE_ENDIAN', &
     iostat     = E_IO)
!
offset = 0
!
! WRITE HEADER
!
write(unit=Unit_VTK,iostat=E_IO)'<?xml version = "1.0"?>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<VTKFile type = "StructuredGrid" version="0.1" byte_order="LittleEndian">'//end_rec
!
write(unit=Unit_VTK,iostat=E_IO)'<StructuredGrid WholeExtent="'  &
                                        //trim(adjustl(cx1))//' '&
                                        //trim(adjustl(cx2))//' '&
                                        //trim(adjustl(cy1))//' '&
                                        //trim(adjustl(cy2))//' '&
                                        //trim(adjustl(cz1))//' '&
                                        //trim(adjustl(cz2))//'">'//end_rec
!
write(unit=Unit_VTK,iostat=E_IO)'<Piece Extent="'           &
                                   //trim(adjustl(cx1))//' '&
                                   //trim(adjustl(cx2))//' '&
                                   //trim(adjustl(cy1))//' '&
                                   //trim(adjustl(cy2))//' '&
                                   //trim(adjustl(cz1))//' '&
                                   //trim(adjustl(cz2))//'">'//end_rec
!
write(unit=Unit_VTK,iostat=E_IO)'<PointData>'//end_rec

  write(coffset,'(I10)')offset
  write(unit=Unit_VTK,iostat=E_IO)'<DataArray type="Float32" Name="kernel beta" NumberOfComponents="1"'&
  &,' format="appended" offset="'//trim(adjustl(coffset))//'" />'//end_rec
  offset = offset+sizeof(KERNEL_values(0,:,:,:))+sizeof(i)
  write(coffset,'(I10)')offset
  write(unit=Unit_VTK,iostat=E_IO)'<DataArray type="Float32" Name="kernel xi" NumberOfComponents="1"'&
  &,' format="appended" offset="'//trim(adjustl(coffset))//'" />'//end_rec
  offset = offset+sizeof(KERNEL_values(0,:,:,:))+sizeof(i)
write(unit=Unit_VTK,iostat=E_IO)'</PointData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<CellData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</CellData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<Points>'//end_rec
write(coffset,'(I10)')offset
write(unit=Unit_VTK,iostat=E_IO)'<DataArray type="Float32" NumberOfComponents="3"'&
&,' format="appended" offset="'//trim(adjustl(coffset))//'" />'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</Points>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</Piece>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</StructuredGrid>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<AppendedData encoding="raw">'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'_'
!
! FREE DISK SPACE FOR DATA 
!

do n = 0,1
  dim_arrays = sizeof(KERNEL_values(n,:,:,:))
  write(unit=Unit_VTK,iostat=E_IO)dim_arrays
  do k = z1,z2
    do j = y1,y2
      do i = x1,x2
        write(unit=Unit_VTK,iostat=E_IO)KERNEL_values(n,i,j,k)
      enddo
    enddo
  enddo
enddo

dim_arrays = sizeof(xc)+sizeof(yc)+sizeof(zc)
write(unit=Unit_VTK,iostat=E_IO)dim_arrays
do k = z1,z2
   do j = y1,y2
      do i = x1,x2
         write(unit=Unit_VTK,iostat=E_IO)xc(i,j,k),yc(i,j,k),zc(i,j,k)
      enddo
   enddo
enddo
write(unit=Unit_VTK,iostat=E_IO)end_rec
!
! WRITE FOOTER
!
write(unit=Unit_VTK,iostat=E_IO)'</AppendedData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</VTKFile>'//end_rec
!
close(unit=Unit_VTK)
  !
  DEALLOCATE(XC,YC,ZC,KERNEL_VALUES)
  !
  endif

endif ! end write vtk

close (31)
if (Tdomain%aniso)   close (32)
call MPI_BARRIER(MPI_COMM_WORLD, code)
if (rg==0) then
    call system("cat kernel_beta_* > kern_beta")
    call system("rm kernel_beta_*")
    if (Tdomain%aniso) then
        call system("cat kernel_xi_* > kern_xi")
        call system("rm kernel_xi_*")
    endif
endif


return
end subroutine save_kernel
