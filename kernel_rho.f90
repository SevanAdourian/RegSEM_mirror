subroutine kernel_rho (Tdomain,dt,ntime)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
doubleprecision, intent(IN) :: dt
integer, intent (IN) :: ntime

integer :: n, x,y,z, ngllx,nglly,ngllz,nf,ne,nv
doubleprecision :: rho
doubleprecision, dimension(:,:,:), allocatable :: kern_rho_loc

do n = 0,Tdomain%n_elem-1
    if (.not. Tdomain%specel(n)%PML) then

ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

allocate (kern_rho_loc(0:ngllx-1,0:nglly-1,0:ngllz-1))

! faces

nf = Tdomain%specel(n)%near_faces(0)
kern_rho_loc(1:ngllx-2,1:nglly-2,0) = Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:nglly-2,0) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:nglly-2,0)&
                                     +Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:nglly-2,1) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:nglly-2,1)&
                                     +Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:nglly-2,2) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:nglly-2,2)
nf = Tdomain%specel(n)%near_faces(1)
kern_rho_loc(1:ngllx-2,0,1:ngllz-2) = Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:ngllz-2,0) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:ngllz-2,0)&
                                     +Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:ngllz-2,1) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:ngllz-2,1)&
                                     +Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:ngllz-2,2) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:ngllz-2,2)
nf = Tdomain%specel(n)%near_faces(2)
kern_rho_loc(ngllx-1,1:nglly-2,1:ngllz-2) = Tdomain%sFace(nf)%sSimu(1)%displ(1:nglly-2,1:ngllz-2,0) * Tdomain%sFace(nf)%sSimu(0)%accel(1:nglly-2,1:ngllz-2,0)&
                                           +Tdomain%sFace(nf)%sSimu(1)%displ(1:nglly-2,1:ngllz-2,1) * Tdomain%sFace(nf)%sSimu(0)%accel(1:nglly-2,1:ngllz-2,1)&
                                           +Tdomain%sFace(nf)%sSimu(1)%displ(1:nglly-2,1:ngllz-2,2) * Tdomain%sFace(nf)%sSimu(0)%accel(1:nglly-2,1:ngllz-2,2)
nf = Tdomain%specel(n)%near_faces(3)
kern_rho_loc(1:ngllx-2,nglly-1,1:ngllz-2) = Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:ngllz-2,0) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:ngllz-2,0)&
                                           +Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:ngllz-2,1) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:ngllz-2,1)&
                                           +Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:ngllz-2,2) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:ngllz-2,2)
nf = Tdomain%specel(n)%near_faces(4)
kern_rho_loc(0,1:nglly-2,1:ngllz-2) = Tdomain%sFace(nf)%sSimu(1)%displ(1:nglly-2,1:ngllz-2,0) * Tdomain%sFace(nf)%sSimu(0)%accel(1:nglly-2,1:ngllz-2,0)&
                                     +Tdomain%sFace(nf)%sSimu(1)%displ(1:nglly-2,1:ngllz-2,1) * Tdomain%sFace(nf)%sSimu(0)%accel(1:nglly-2,1:ngllz-2,1)&
                                     +Tdomain%sFace(nf)%sSimu(1)%displ(1:nglly-2,1:ngllz-2,2) * Tdomain%sFace(nf)%sSimu(0)%accel(1:nglly-2,1:ngllz-2,2)
nf = Tdomain%specel(n)%near_faces(5)
kern_rho_loc(1:ngllx-2,1:nglly-2,ngllz-1) = Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:nglly-2,0) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:nglly-2,0)&
                                           +Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:nglly-2,1) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:nglly-2,1)&
                                           +Tdomain%sFace(nf)%sSimu(1)%displ(1:ngllx-2,1:nglly-2,2) * Tdomain%sFace(nf)%sSimu(0)%accel(1:ngllx-2,1:nglly-2,2)

! edges


ne = Tdomain%specel(n)%near_edges(0)
kern_rho_loc(1:ngllx-2,0,0) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,0)&
                             +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,1)&
                             +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,2)
ne = Tdomain%specel(n)%near_edges(1)
kern_rho_loc(ngllx-1,1:nglly-2,0) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,0)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,1)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,2)
ne = Tdomain%specel(n)%near_edges(2)
kern_rho_loc(1:ngllx-2,nglly-1,0) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,0)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,1)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,2)
ne = Tdomain%specel(n)%near_edges(3)
kern_rho_loc(0,1:nglly-2,0) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,0)&
                             +Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,1)&
                             +Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,2)
ne = Tdomain%specel(n)%near_edges(4)
kern_rho_loc(ngllx-1,0,1:ngllz-2) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,0)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,1)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,2)
ne = Tdomain%specel(n)%near_edges(5)
kern_rho_loc(1:ngllx-2,0,ngllz-1) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,0)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,1)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,2)
ne = Tdomain%specel(n)%near_edges(6)
kern_rho_loc(0,0,1:ngllz-2) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,0)&
                             +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,1)&
                             +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,2)
ne = Tdomain%specel(n)%near_edges(7)
kern_rho_loc(ngllx-1,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,0)&
                                         +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,1)&
                                         +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,2)
ne = Tdomain%specel(n)%near_edges(8)
kern_rho_loc(ngllx-1,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,0)&
                                         +Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,1)&
                                         +Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,2)
ne = Tdomain%specel(n)%near_edges(9)
kern_rho_loc(1:ngllx-2,nglly-1,ngllz-1) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,0)&
                                         +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,1)&
                                         +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllx-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllx-2,2)
ne = Tdomain%specel(n)%near_edges(10)
kern_rho_loc(0,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,0)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,1)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:ngllz-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:ngllz-2,2)
ne = Tdomain%specel(n)%near_edges(11)
kern_rho_loc(0,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,0)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,0)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,1)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,1)&
                                   +Tdomain%sEdge(ne)%sSimu(1)%displ(1:nglly-2,2)*Tdomain%sEdge(ne)%sSimu(0)%accel(1:nglly-2,2)

! vertices


nv = Tdomain%specel(n)%near_vertices(0)
kern_rho_loc(0,0,0) = Tdomain%sVertex(nv)%sSimu(1)%displ(0)*Tdomain%sVertex(nv)%sSimu(0)%accel(0)& 
                     +Tdomain%sVertex(nv)%sSimu(1)%displ(1)*Tdomain%sVertex(nv)%sSimu(0)%accel(1)& 
                     +Tdomain%sVertex(nv)%sSimu(1)%displ(2)*Tdomain%sVertex(nv)%sSimu(0)%accel(2)
nv = Tdomain%specel(n)%near_vertices(1)
kern_rho_loc(ngllx-1,0,0) = Tdomain%sVertex(nv)%sSimu(1)%displ(0)*Tdomain%sVertex(nv)%sSimu(0)%accel(0)& 
                           +Tdomain%sVertex(nv)%sSimu(1)%displ(1)*Tdomain%sVertex(nv)%sSimu(0)%accel(1)& 
                           +Tdomain%sVertex(nv)%sSimu(1)%displ(2)*Tdomain%sVertex(nv)%sSimu(0)%accel(2)
nv = Tdomain%specel(n)%near_vertices(2)
kern_rho_loc(ngllx-1,nglly-1,0) = Tdomain%sVertex(nv)%sSimu(1)%displ(0)*Tdomain%sVertex(nv)%sSimu(0)%accel(0)& 
                                 +Tdomain%sVertex(nv)%sSimu(1)%displ(1)*Tdomain%sVertex(nv)%sSimu(0)%accel(1)& 
                                 +Tdomain%sVertex(nv)%sSimu(1)%displ(2)*Tdomain%sVertex(nv)%sSimu(0)%accel(2)
nv = Tdomain%specel(n)%near_vertices(3)
kern_rho_loc(0,nglly-1,0) = Tdomain%sVertex(nv)%sSimu(1)%displ(0)*Tdomain%sVertex(nv)%sSimu(0)%accel(0)& 
                           +Tdomain%sVertex(nv)%sSimu(1)%displ(1)*Tdomain%sVertex(nv)%sSimu(0)%accel(1)& 
                           +Tdomain%sVertex(nv)%sSimu(1)%displ(2)*Tdomain%sVertex(nv)%sSimu(0)%accel(2)
nv = Tdomain%specel(n)%near_vertices(4)
kern_rho_loc(0,0,ngllz-1) = Tdomain%sVertex(nv)%sSimu(1)%displ(0)*Tdomain%sVertex(nv)%sSimu(0)%accel(0)& 
                           +Tdomain%sVertex(nv)%sSimu(1)%displ(1)*Tdomain%sVertex(nv)%sSimu(0)%accel(1)& 
                           +Tdomain%sVertex(nv)%sSimu(1)%displ(2)*Tdomain%sVertex(nv)%sSimu(0)%accel(2)
nv = Tdomain%specel(n)%near_vertices(5)
kern_rho_loc(ngllx-1,0,ngllz-1) = Tdomain%sVertex(nv)%sSimu(1)%displ(0)*Tdomain%sVertex(nv)%sSimu(0)%accel(0)& 
                                 +Tdomain%sVertex(nv)%sSimu(1)%displ(1)*Tdomain%sVertex(nv)%sSimu(0)%accel(1)& 
                                 +Tdomain%sVertex(nv)%sSimu(1)%displ(2)*Tdomain%sVertex(nv)%sSimu(0)%accel(2)
nv = Tdomain%specel(n)%near_vertices(6)
kern_rho_loc(ngllx-1,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%sSimu(1)%displ(0)*Tdomain%sVertex(nv)%sSimu(0)%accel(0)& 
                                       +Tdomain%sVertex(nv)%sSimu(1)%displ(1)*Tdomain%sVertex(nv)%sSimu(0)%accel(1)& 
                                       +Tdomain%sVertex(nv)%sSimu(1)%displ(2)*Tdomain%sVertex(nv)%sSimu(0)%accel(2)
nv = Tdomain%specel(n)%near_vertices(7)
kern_rho_loc(0,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%sSimu(1)%displ(0)*Tdomain%sVertex(nv)%sSimu(0)%accel(0)& 
                                 +Tdomain%sVertex(nv)%sSimu(1)%displ(1)*Tdomain%sVertex(nv)%sSimu(0)%accel(1)& 
                                 +Tdomain%sVertex(nv)%sSimu(1)%displ(2)*Tdomain%sVertex(nv)%sSimu(0)%accel(2)

! interior nodes


kern_rho_loc(1:ngllx-2,1:nglly-2,1:ngllz-2) =  Tdomain%specel(n)%sSimu(1)%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0)*Tdomain%specel(n)%sSimu(0)%Accel(1:ngllx-2,1:nglly-2,1:ngllz-2,0) &
                                             + Tdomain%specel(n)%sSimu(1)%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,1)*Tdomain%specel(n)%sSimu(0)%Accel(1:ngllx-2,1:nglly-2,1:ngllz-2,1) &
                                             + Tdomain%specel(n)%sSimu(1)%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,2)*Tdomain%specel(n)%sSimu(0)%Accel(1:ngllx-2,1:nglly-2,1:ngllz-2,2) 
! multiply by density

kern_rho_loc(:,:,:) = kern_rho_loc(:,:,:) * Tdomain%specel(n)%Density(:,:,:)
Tdomain%specel(n)%Kern_rho(:,:,:) = Tdomain%specel(n)%Kern_rho(:,:,:) + dt*kern_rho_loc(:,:,:)
!write(*,*) "kern rho is", kern_rho_loc(:,:,:) 
! old version by paul cupillard
!        do z = 1,ngllz-2
!         do y = 1,nglly-2
!          do x = 1,ngllx-2
!              rho = Tdomain%specel(n)%Density(x,y,z)
!              kern_rho_loc(x,y,z) = -rho * &
!                                    ( Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,0)*Tdomain%specel(n)%sSimu(0)%Accel(x,y,z,0) &
!                                    + Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,1)*Tdomain%specel(n)%sSimu(0)%Accel(x,y,z,1) &
!                                    + Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,2)*Tdomain%specel(n)%sSimu(0)%Accel(x,y,z,2) )
!              Tdomain%specel(n)%Kern_rho(x,y,z) = Tdomain%specel(n)%Kern_rho(x,y,z) + dt*kern_rho_loc(x,y,z)
!          enddo
!         enddo
!        enddo
        !!! Reste a chopper les valeurs des differents champs sur les faces, les edges et les vertices !!!
        deallocate (kern_rho_loc)
    endif
enddo


return
end subroutine kernel_rho
