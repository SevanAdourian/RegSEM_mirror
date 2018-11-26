!=================================================
subroutine impose_init (Tdomain, ntime, ninit, rg)
!=================================================
!
! Written by Yder MASSON (masson@ipgp.fr) june 20 2012  
! 
! reset fields variables (point matching) when running the time reversal mirror simulation with attenation 
!
use sdomains
!
implicit none
!
type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: ntime, ninit, rg
!
integer :: len_DP, n, ngllx, nglly, ngllz, comp, nf, ne, nv, n_solid, i, ngll, ngll1, ngll2, mypos, recl_init
logical :: pml
character*60 :: init_file
real, dimension (:,:,:,:), allocatable :: win_mirror_loc
logical, dimension(:), allocatable :: L_Face, L_Edge, L_Vertex
!
if (ntime==0) then
!
len_DP = 0
!
do n = 0,Tdomain%n_elem-1
    if (.not.Tdomain%specel(n)%PML) then
        len_DP = len_DP + sizeof(Tdomain%specel(n)%sSimu(0)%Displ)+sizeof(Tdomain%specel(n)%sSimu(0)%Veloc)
        if(Tdomain%n_sls>0)then
        len_DP = len_DP                                   &
        &+sizeof(Tdomain%specel(n)%sSimu(0)%R_xx_)         &
        &+sizeof(Tdomain%specel(n)%sSimu(0)%R_yy_)         &
        &+sizeof(Tdomain%specel(n)%sSimu(0)%R_xy_)         &
        &+sizeof(Tdomain%specel(n)%sSimu(0)%R_xz_)         &
        &+sizeof(Tdomain%specel(n)%sSimu(0)%R_yz_)         & 
        &+sizeof(Tdomain%specel(n)%sSimu(0)%epsilondev_xx_)&
        &+sizeof(Tdomain%specel(n)%sSimu(0)%epsilondev_yy_)&
        &+sizeof(Tdomain%specel(n)%sSimu(0)%epsilondev_xy_)&
        &+sizeof(Tdomain%specel(n)%sSimu(0)%epsilondev_xz_)&
        &+sizeof(Tdomain%specel(n)%sSimu(0)%epsilondev_yz_)
        endif
    endif
enddo
!
do n = 0,Tdomain%n_face-1
    if (.not.Tdomain%sFace(n)%PML) then 
        len_DP = len_DP + sizeof(Tdomain%sface(n)%sSimu(0)%Displ)+sizeof(Tdomain%sface(n)%sSimu(0)%Veloc)
    endif
enddo
!
do n = 0,Tdomain%n_edge-1
    if (.not.Tdomain%sEdge(n)%PML) then
        len_DP = len_DP + sizeof(Tdomain%sedge(n)%sSimu(0)%Displ)+sizeof(Tdomain%sedge(n)%sSimu(0)%Veloc)
endif
!
enddo
do n = 0,Tdomain%n_vertex-1
    if (.not.Tdomain%sVertex(n)%PML) then
        len_DP = len_DP + sizeof(Tdomain%svertex(n)%sSimu(0)%Displ)+sizeof(Tdomain%svertex(n)%sSimu(0)%Veloc)
    endif
enddo
!
Tdomain%recl_init = len_DP+sizeof(mypos)
!   
   if (len_DP/=0) then
      write (init_file,"(a,I3.3)") "init",rg
      open (45,file=init_file,access='stream',form="unformatted",status='old')
   endif
!
endif
!
if(mod(ntime,ninit)==0)then
!
! Zero fields
!
n_solid = Tdomain%n_sls
!
do n = 0,Tdomain%n_elem-1
    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz
    Tdomain%specel(n)%ssimu(0)%Veloc = 0.d0
    Tdomain%specel(n)%ssimu(0)%Accel = 0.d0
    if (Tdomain%specel(n)%PML) then
            Tdomain%specel(n)%ssimu(0)%Diagonal_Stress = 0.d0
            Tdomain%specel(n)%ssimu(0)%Diagonal_Stress1 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Diagonal_Stress2 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Diagonal_Stress3 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Residual_Stress = 0.d0
            Tdomain%specel(n)%ssimu(0)%Residual_Stress1 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Residual_Stress2 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Veloc1 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Veloc2 = 0.d0
            Tdomain%specel(n)%ssimu(0)%Veloc3 = 0.d0
        if (Tdomain%curve) then
                Tdomain%specel(n)%ssimu(0)%Residual_Stress3 = 0.d0
        endif
    else
            Tdomain%specel(n)%ssimu(0)%Displ = 0.d0
            Tdomain%specel(n)%ssimu(0)%V0 = 0.d0
        if (n_solid>0) then
                Tdomain%specel(n)%ssimu(0)%epsilondev_xx_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%epsilondev_yy_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%epsilondev_xy_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%epsilondev_xz_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%epsilondev_yz_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_xx_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_yy_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_xy_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_xz_ = 0.d0
                Tdomain%specel(n)%ssimu(0)%R_yz_ = 0.d0
        endif
    endif
enddo

do n = 0,Tdomain%n_face-1
    ngll1 = Tdomain%sFace(n)%ngll1
    ngll2 = Tdomain%sFace(n)%ngll2
        Tdomain%sFace(n)%ssimu(0)%Veloc = 0.d0
        Tdomain%sFace(n)%ssimu(0)%Accel = 0.d0
        Tdomain%sFace(n)%ssimu(0)%Forces = 0.d0
    if (Tdomain%sFace(n)%PML) then 
            Tdomain%sFace(n)%ssimu(0)%Veloc1 = 0.d0
            Tdomain%sFace(n)%ssimu(0)%Veloc2 = 0.d0
            Tdomain%sFace(n)%ssimu(0)%Veloc3 = 0.d0
    else
            Tdomain%sFace(n)%ssimu(0)%Displ = 0.d0
            Tdomain%sFace(n)%ssimu(0)%V0 = 0.d0
    endif
enddo

do n = 0,Tdomain%n_edge-1
    ngll = Tdomain%sEdge(n)%ngll
        Tdomain%sEdge(n)%ssimu(0)%Veloc = 0.d0
        Tdomain%sEdge(n)%ssimu(0)%Accel = 0.d0
        Tdomain%sEdge(n)%ssimu(0)%Forces = 0.d0
    if (Tdomain%sEdge(n)%PML) then
            Tdomain%sEdge(n)%ssimu(0)%Veloc1 = 0.d0
            Tdomain%sEdge(n)%ssimu(0)%Veloc2 = 0.d0
            Tdomain%sEdge(n)%ssimu(0)%Veloc3 = 0.d0
    else
            Tdomain%sEdge(n)%ssimu(0)%Displ = 0.d0
            Tdomain%sEdge(n)%ssimu(0)%V0 = 0.d0
    endif
enddo

do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%ssimu(0)%Veloc = 0.d0
        Tdomain%sVertex(n)%ssimu(0)%Accel = 0.d0
        Tdomain%sVertex(n)%ssimu(0)%Forces = 0.d0
    if (Tdomain%sVertex(n)%PML) then
            Tdomain%sVertex(n)%ssimu(0)%Veloc1 = 0.d0
            Tdomain%sVertex(n)%ssimu(0)%Veloc2 = 0.d0
            Tdomain%sVertex(n)%ssimu(0)%Veloc3 = 0.d0
    else
            Tdomain%sVertex(n)%ssimu(0)%Displ = 0.d0
            Tdomain%sVertex(n)%ssimu(0)%V0 = 0.d0
    endif
enddo
!
mypos = ntime/ninit*Tdomain%recl_init+1 
read(45,pos=mypos)recl_init
if(Tdomain%recl_init/=recl_init)write(*,*)"Warning!!! record length dosen't match in rec_init/impose_init"
!
do n = 0,Tdomain%n_elem-1
    if (.not.Tdomain%specel(n)%PML) then
        read (45)Tdomain%specel(n)%sSimu(0)%Displ,Tdomain%specel(n)%sSimu(0)%Veloc
        if(Tdomain%n_sls>0)then
        read (45)                                 &
        &Tdomain%specel(n)%sSimu(0)%R_xx_,         &
        &Tdomain%specel(n)%sSimu(0)%R_yy_,         &
        &Tdomain%specel(n)%sSimu(0)%R_xy_,         &
        &Tdomain%specel(n)%sSimu(0)%R_xz_,         &
        &Tdomain%specel(n)%sSimu(0)%R_yz_,         & 
        &Tdomain%specel(n)%sSimu(0)%epsilondev_xx_,&
        &Tdomain%specel(n)%sSimu(0)%epsilondev_yy_,&
        &Tdomain%specel(n)%sSimu(0)%epsilondev_xy_,&
        &Tdomain%specel(n)%sSimu(0)%epsilondev_xz_,&
        &Tdomain%specel(n)%sSimu(0)%epsilondev_yz_
        endif
    endif
enddo
!
do n = 0,Tdomain%n_face-1
    if (.not.Tdomain%sFace(n)%PML) then 
        read (45)Tdomain%sface(n)%sSimu(0)%Displ,Tdomain%sface(n)%sSimu(0)%Veloc
    endif
enddo
!
do n = 0,Tdomain%n_edge-1
    if (.not.Tdomain%sEdge(n)%PML) then
        read (45)Tdomain%sedge(n)%sSimu(0)%Displ,Tdomain%sedge(n)%sSimu(0)%Veloc
    endif
enddo
!
do n = 0,Tdomain%n_vertex-1
    if (.not.Tdomain%sVertex(n)%PML) then
        read (45)Tdomain%svertex(n)%sSimu(0)%Displ,Tdomain%svertex(n)%sSimu(0)%Veloc
    endif
enddo
!
!   read (45,rec=ntime/ninit+1) &
!        &(Tdomain%specel(n)%sSimu(0)%Displ,n=0,Tdomain%n_elem-1),         &
!        &(Tdomain%sFace(n)%sSimu(0)%Displ,n=0,Tdomain%n_face-1),          &
!        &(Tdomain%sEdge(n)%sSimu(0)%Displ,n=0,Tdomain%n_edge-1),          &
!        &(Tdomain%sVertex(n)%sSimu(0)%Displ,n=0,Tdomain%n_vertex-1),      &
!        &                                                                 &
!        &(Tdomain%specel(n)%sSimu(0)%Veloc,n=0,Tdomain%n_elem-1),         &
!        &(Tdomain%sFace(n)%sSimu(0)%Veloc,n=0,Tdomain%n_face-1),          &
!        &(Tdomain%sEdge(n)%sSimu(0)%Veloc,n=0,Tdomain%n_edge-1),          &
!        &(Tdomain%sVertex(n)%sSimu(0)%Veloc,n=0,Tdomain%n_vertex-1),      &
!        &                                                                 &
!        &(Tdomain%specel(n)%sSimu(0)%R_xx_,n=0,Tdomain%n_elem-1),         &
!        &(Tdomain%specel(n)%sSimu(0)%R_yy_,n=0,Tdomain%n_elem-1),         &
!        &(Tdomain%specel(n)%sSimu(0)%R_xy_,n=0,Tdomain%n_elem-1),         &
!        &(Tdomain%specel(n)%sSimu(0)%R_xz_,n=0,Tdomain%n_elem-1),         &
!        &(Tdomain%specel(n)%sSimu(0)%R_yz_,n=0,Tdomain%n_elem-1),         &
!        &                                                                 &
!        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xx_,n=0,Tdomain%n_elem-1),&
!        &(Tdomain%specel(n)%sSimu(0)%epsilondev_yy_,n=0,Tdomain%n_elem-1),&
!        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xy_,n=0,Tdomain%n_elem-1),&
!        &(Tdomain%specel(n)%sSimu(0)%epsilondev_xz_,n=0,Tdomain%n_elem-1),&
!        &(Tdomain%specel(n)%sSimu(0)%epsilondev_yz_,n=0,Tdomain%n_elem-1)
!
! Apply mirror to initial conditions
!
 allocate (L_Face(0:Tdomain%n_face-1))
 L_Face = .true.
 allocate (L_Edge(0:Tdomain%n_edge-1))
 L_Edge = .true.
 allocate (L_Vertex(0:Tdomain%n_vertex-1))
 L_Vertex = .true.

do n = 0,Tdomain%n_elem-1

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    allocate (win_mirror_loc(0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:2))

    if(Tdomain%specel(n)%mirror_position == 0)then
        win_mirror_loc(:,:,:,:) = 0
    elseif(Tdomain%specel(n)%mirror_position == 2)then
        win_mirror_loc(:,:,:,:) = 1
    else
        win_mirror_loc(:,:,:,0) = Tdomain%specel(n)%win_mirror(:,:,:)
        win_mirror_loc(:,:,:,1) = Tdomain%specel(n)%win_mirror(:,:,:)
        win_mirror_loc(:,:,:,2) = Tdomain%specel(n)%win_mirror(:,:,:)
    endif

    if (.not.Tdomain%specel(n)%PML) then

    Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,0) = Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,0) * win_mirror_loc(1:ngllx-2, 1:nglly-2, 1:ngllz-2,0)
    Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,1) = Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,1) * win_mirror_loc(1:ngllx-2, 1:nglly-2, 1:ngllz-2,1)
    Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,2) = Tdomain%specel(n)%sSimu(0)%Displ(:,:,:,2) * win_mirror_loc(1:ngllx-2, 1:nglly-2, 1:ngllz-2,2)

    Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,0) = Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,0) * win_mirror_loc(1:ngllx-2, 1:nglly-2, 1:ngllz-2,0)
    Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,1) = Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,1) * win_mirror_loc(1:ngllx-2, 1:nglly-2, 1:ngllz-2,1)
    Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,2) = Tdomain%specel(n)%sSimu(0)%Veloc(:,:,:,2) * win_mirror_loc(1:ngllx-2, 1:nglly-2, 1:ngllz-2,2)

    endif
!
nf = Tdomain%specel(n)%near_faces(0)
pml = .not.Tdomain%sFace(nf)%PML
if(pml.and.L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,:) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,:)&
                                                                  *win_mirror_loc(1:ngllx-2,1:nglly-2,0,:)

nf = Tdomain%specel(n)%near_faces(1)
pml = .not.Tdomain%sFace(nf)%PML
if(pml.and.L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:ngllz-2,:) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:ngllz-2,:)&
                                                                  *win_mirror_loc(1:ngllx-2,0,1:ngllz-2,:)

nf = Tdomain%specel(n)%near_faces(2)
pml = .not.Tdomain%sFace(nf)%PML
if(pml.and.L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Displ(1:nglly-2,1:ngllz-2,:) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:nglly-2,1:ngllz-2,:)&
                                                                  *win_mirror_loc(ngllx-1,1:nglly-2,1:ngllz-2,:)

nf = Tdomain%specel(n)%near_faces(3)
pml = .not.Tdomain%sFace(nf)%PML
if(pml.and.L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:ngllz-2,:) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:ngllz-2,:)&
                                                                  *win_mirror_loc(1:ngllx-2,nglly-1,1:ngllz-2,:)

nf = Tdomain%specel(n)%near_faces(4)
pml = .not.Tdomain%sFace(nf)%PML
if(pml.and.L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Displ(1:nglly-2,1:ngllz-2,:) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:nglly-2,1:ngllz-2,:)&
                                                                  *win_mirror_loc(0,1:nglly-2,1:ngllz-2,:)

nf = Tdomain%specel(n)%near_faces(5)
pml = .not.Tdomain%sFace(nf)%PML
if(pml.and.L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,:) = Tdomain%sFace(nf)%sSimu(0)%Displ(1:ngllx-2,1:nglly-2,:)&
                                                                  *win_mirror_loc(1:ngllx-2,1:nglly-2,ngllz-1,:)
!

ne = Tdomain%specel(n)%near_edges(0)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,:)&
                                                        *win_mirror_loc(1:ngllx-2,0,0,:) 

ne = Tdomain%specel(n)%near_edges(1)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,:)&
                                                         *win_mirror_loc(ngllx-1,1:nglly-2,0,:)

ne = Tdomain%specel(n)%near_edges(2)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,:)&
                                                        *win_mirror_loc(1:ngllx-2,nglly-1,0,:) 

ne = Tdomain%specel(n)%near_edges(3)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,:)&
                                                         *win_mirror_loc(0,1:nglly-2,0,:) 

ne = Tdomain%specel(n)%near_edges(4)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,:)&
                                                        *win_mirror_loc(ngllx-1,0,1:ngllz-2,:) 

ne = Tdomain%specel(n)%near_edges(5)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,:)&
                                                        *win_mirror_loc(1:ngllx-2,0,ngllz-1,:) 

ne = Tdomain%specel(n)%near_edges(6)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,:)&
                                                        *win_mirror_loc(0,0,1:ngllz-2,:) 

ne = Tdomain%specel(n)%near_edges(7)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,:)&
                                                         * win_mirror_loc(ngllx-1,nglly-1,1:ngllz-2,:)

ne = Tdomain%specel(n)%near_edges(8)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,:)&
                                                         *win_mirror_loc(ngllx-1,1:nglly-2,ngllz-1,:)

ne = Tdomain%specel(n)%near_edges(9)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllx-2,:)&
                                                        *win_mirror_loc(1:ngllx-2,nglly-1,ngllz-1,:) 

ne = Tdomain%specel(n)%near_edges(10)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:ngllz-2,:)&
                                                        *win_mirror_loc(0,nglly-1,1:ngllz-2,:) 

ne = Tdomain%specel(n)%near_edges(11)
pml = .not.Tdomain%sEdge(ne)%PML
if(pml.and.L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Displ(1:nglly-2,:)&
                                                        *win_mirror_loc(0,1:nglly-2,ngllz-1,:) 
!

nv = Tdomain%specel(n)%near_vertices(0)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml.and.L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Displ(:)= Tdomain%sVertex(nv)%sSimu(0)%Displ(:)&
                                                *win_mirror_loc(0,0,0,:) 

nv = Tdomain%specel(n)%near_vertices(1)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml.and.L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Displ(:) = Tdomain%sVertex(nv)%sSimu(0)%Displ(:)&
                                                *win_mirror_loc(ngllx-1,0,0,:)

nv = Tdomain%specel(n)%near_vertices(2)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml.and.L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Displ(:) = Tdomain%sVertex(nv)%sSimu(0)%Displ(:)&
                                                *win_mirror_loc(ngllx-1,nglly-1,0,:) 

nv = Tdomain%specel(n)%near_vertices(3)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml.and.L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Displ(:) = Tdomain%sVertex(nv)%sSimu(0)%Displ(:)&
                                                *win_mirror_loc(0,nglly-1,0,:)

nv = Tdomain%specel(n)%near_vertices(4)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml.and.L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Displ(:) = Tdomain%sVertex(nv)%sSimu(0)%Displ(:)&
                                                 *win_mirror_loc(0,0,ngllz-1,:)

nv = Tdomain%specel(n)%near_vertices(5)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml.and.L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Displ(:) = Tdomain%sVertex(nv)%sSimu(0)%Displ(:)&
                                                 *win_mirror_loc(ngllx-1,0,ngllz-1,:)

nv = Tdomain%specel(n)%near_vertices(6)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml.and.L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Displ(:) = Tdomain%sVertex(nv)%sSimu(0)%Displ(:)&
                                                 *win_mirror_loc(ngllx-1,nglly-1,ngllz-1,:)

nv = Tdomain%specel(n)%near_vertices(7)
pml = .not.Tdomain%sVertex(nv)%PML
if(pml.and.L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Displ(:) = Tdomain%sVertex(nv)%sSimu(0)%Displ(:)&
                                                *win_mirror_loc(0,nglly-1,ngllz-1,:)
!
! veloc
!
nf = Tdomain%specel(n)%near_faces(0)
if(L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,:) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,:)&
                                                                  *win_mirror_loc(1:ngllx-2,1:nglly-2,0,:)
L_Face(nf) = .false.

nf = Tdomain%specel(n)%near_faces(1)
if(L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:ngllz-2,:) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:ngllz-2,:)&
                                                                  *win_mirror_loc(1:ngllx-2,0,1:ngllz-2,:)
L_Face(nf) = .false.

nf = Tdomain%specel(n)%near_faces(2)
if(L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Veloc(1:nglly-2,1:ngllz-2,:) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:nglly-2,1:ngllz-2,:)&
                                                                  *win_mirror_loc(ngllx-1,1:nglly-2,1:ngllz-2,:)
L_Face(nf) = .false.

nf = Tdomain%specel(n)%near_faces(3)
if(L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:ngllz-2,:) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:ngllz-2,:)&
                                                                  *win_mirror_loc(1:ngllx-2,nglly-1,1:ngllz-2,:)
L_Face(nf) = .false.

nf = Tdomain%specel(n)%near_faces(4)
if(L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Veloc(1:nglly-2,1:ngllz-2,:) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:nglly-2,1:ngllz-2,:)&
                                                                  *win_mirror_loc(0,1:nglly-2,1:ngllz-2,:)
L_Face(nf) = .false.

nf = Tdomain%specel(n)%near_faces(5)
if(L_Face(nf))Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,:) = Tdomain%sFace(nf)%sSimu(0)%Veloc(1:ngllx-2,1:nglly-2,:)&
                                                                  *win_mirror_loc(1:ngllx-2,1:nglly-2,ngllz-1,:)
L_Face(nf) = .false.
                                                                  !

ne = Tdomain%specel(n)%near_edges(0)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,:)&
                                                        *win_mirror_loc(1:ngllx-2,0,0,:) 
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(1)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,:)&
                                                         *win_mirror_loc(ngllx-1,1:nglly-2,0,:)
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(2)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,:)&
                                                        *win_mirror_loc(1:ngllx-2,nglly-1,0,:) 
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(3)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,:)&
                                                         *win_mirror_loc(0,1:nglly-2,0,:) 
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(4)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,:)&
                                                        *win_mirror_loc(ngllx-1,0,1:ngllz-2,:) 
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(5)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,:)&
                                                        *win_mirror_loc(1:ngllx-2,0,ngllz-1,:) 
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(6)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,:)&
                                                        *win_mirror_loc(0,0,1:ngllz-2,:) 
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(7)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,:)&
                                                         * win_mirror_loc(ngllx-1,nglly-1,1:ngllz-2,:)
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(8)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,:)&
                                                         *win_mirror_loc(ngllx-1,1:nglly-2,ngllz-1,:)
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(9)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllx-2,:)&
                                                        *win_mirror_loc(1:ngllx-2,nglly-1,ngllz-1,:) 
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(10)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:ngllz-2,:)&
                                                        *win_mirror_loc(0,nglly-1,1:ngllz-2,:) 
L_Edge(ne) = .false.

ne = Tdomain%specel(n)%near_edges(11)
if(L_Edge(ne))Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,:) = Tdomain%sEdge(ne)%sSimu(0)%Veloc(1:nglly-2,:)&
                                                        *win_mirror_loc(0,1:nglly-2,ngllz-1,:) 
L_Edge(ne) = .false.
!

nv = Tdomain%specel(n)%near_vertices(0)
if(L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)= Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)&
                                                *win_mirror_loc(0,0,0,:) 
L_Vertex(nv) = .false.

nv = Tdomain%specel(n)%near_vertices(1)
if(L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Veloc(:) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)&
                                                *win_mirror_loc(ngllx-1,0,0,:)
L_Vertex(nv) = .false.

nv = Tdomain%specel(n)%near_vertices(2)
if(L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Veloc(:) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)&
                                                *win_mirror_loc(ngllx-1,nglly-1,0,:) 
L_Vertex(nv) = .false.

nv = Tdomain%specel(n)%near_vertices(3)
if(L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Veloc(:) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)&
                                                *win_mirror_loc(0,nglly-1,0,:)
L_Vertex(nv) = .false.

nv = Tdomain%specel(n)%near_vertices(4)
if(L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Veloc(:) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)&
                                                 *win_mirror_loc(0,0,ngllz-1,:)
L_Vertex(nv) = .false.

nv = Tdomain%specel(n)%near_vertices(5)
if(L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Veloc(:) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)&
                                                 *win_mirror_loc(ngllx-1,0,ngllz-1,:)
L_Vertex(nv) = .false.

nv = Tdomain%specel(n)%near_vertices(6)
if(L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Veloc(:) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)&
                                                 *win_mirror_loc(ngllx-1,nglly-1,ngllz-1,:)
L_Vertex(nv) = .false.

nv = Tdomain%specel(n)%near_vertices(7)
if(L_Vertex(nv))Tdomain%sVertex(nv)%sSimu(0)%Veloc(:) = Tdomain%sVertex(nv)%sSimu(0)%Veloc(:)&
                                                *win_mirror_loc(0,nglly-1,ngllz-1,:)
L_Vertex(nv) = .false.
!
deallocate (win_mirror_loc)
!
enddo
!
deallocate(L_Face,L_Vertex,L_Edge)
!
endif
!
if (ntime==Tdomain%sTimeParam%ntime-1)   close (45)
!
!=========================
end subroutine impose_init
!=========================
