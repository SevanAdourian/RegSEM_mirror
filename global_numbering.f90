subroutine global_numbering (Tdomain)

use sdomains

implicit none

include 'mpif.h'

type(domain), target, intent (INOUT) :: Tdomain

integer :: n, icount, i,j,k, ngllx,nglly,ngllz, nf,ne,nv,code
integer, dimension(:,:), allocatable :: indx_gll_all,indy_gll_all,indz_gll_all
!
! Global numbering common to all proc
!
! we compute :
!
! the global indexes (i.e. in the entire gll structured mesh) of the first gll in each element
! if n is the local element number, we have : 
! indx_gll(Tdomain%specel(n)%indx)
! indy_gll(Tdomain%specel(n)%indy)
! indy_gll(Tdomain%specel(n)%indz)
! that are the 3D indexes of the first gll
! where
! Tdomain%specel(n)%indx 
! Tdomain%specel(n)%indy  
! Tdomain%specel(n)%indz 
! are the 3D indexes in the entire structured element mesh
!
! the dimension : ngll_gridx, ngll_grid_y and ngll_grid_z
! are the total numer of glls in the 3 spatial dimension
!
! using the indexes above, it is possible toi assign each gll a unique 1D (global number) 
! index as follow :
!
! say we have
!
! i = indx_gll(Tdomain%specel(n)%indx)
! j = indy_gll(Tdomain%specel(n)%indy) 
! k = indz_gll(Tdomain%specel(n)%indz)
! (i.e the indexes of the first corner gll in the current element)
!
! and
!
! igll, jgll and kgll are the indexes of the gll node in the current element
!
! and
!
! nx = ngll_grid_x
! ny = ngll_grid_y
! nz = ngll_grid_z
!
! then (we assume i,j,k, igll, jgll and kgll starts at zero)
!
! 1D_index = (i+igll)+(j+jgll)*nx+(k+kgll)*(nx*ny)
!
allocate(Tdomain%indx_gll(0:Tdomain%n_elem_x-1))
allocate(Tdomain%indy_gll(0:Tdomain%n_elem_y-1))
allocate(Tdomain%indz_gll(0:Tdomain%n_elem_z-1))
!
Tdomain%indx_gll = 0
Tdomain%indy_gll = 0
Tdomain%indz_gll = 0
!
do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   Tdomain%indx_gll(Tdomain%specel(n)%indx) = ngllx-1
   Tdomain%indy_gll(Tdomain%specel(n)%indy) = nglly-1
   Tdomain%indz_gll(Tdomain%specel(n)%indz) = ngllz-1
enddo
!
allocate(indx_gll_all(0:Tdomain%n_elem_x-1,Tdomain%n_proc))
allocate(indy_gll_all(0:Tdomain%n_elem_y-1,Tdomain%n_proc))
allocate(indz_gll_all(0:Tdomain%n_elem_z-1,Tdomain%n_proc))
!
call mpi_allgather(Tdomain%indx_gll(0),Tdomain%n_elem_x,mpi_integer,indx_gll_all(0,1),Tdomain%n_elem_x,mpi_integer,mpi_comm_world,code)
call mpi_allgather(Tdomain%indy_gll(0),Tdomain%n_elem_y,mpi_integer,indy_gll_all(0,1),Tdomain%n_elem_y,mpi_integer,mpi_comm_world,code)
call mpi_allgather(Tdomain%indz_gll(0),Tdomain%n_elem_z,mpi_integer,indz_gll_all(0,1),Tdomain%n_elem_z,mpi_integer,mpi_comm_world,code)
!
do n = 0,Tdomain%n_elem_x-1
   Tdomain%indx_gll(n) = maxval(indx_gll_all(n,:))        
enddo
do n = 0,Tdomain%n_elem_y-1
   Tdomain%indy_gll(n) = maxval(indy_gll_all(n,:))        
enddo
do n = 0,Tdomain%n_elem_z-1
   Tdomain%indz_gll(n) = maxval(indz_gll_all(n,:))        
enddo
!
Tdomain%ngll_grid_x = sum(Tdomain%indx_gll)+1
Tdomain%ngll_grid_y = sum(Tdomain%indy_gll)+1 
Tdomain%ngll_grid_z = sum(Tdomain%indz_gll)+1 
!
! cumulative gll count in the 3 dimension
!
do n = 1,Tdomain%n_elem_x-1
   Tdomain%indx_gll(n) = Tdomain%indx_gll(n)+Tdomain%indx_gll(n-1)        
enddo
do n = 1,Tdomain%n_elem_y-1
   Tdomain%indy_gll(n) = Tdomain%indy_gll(n)+Tdomain%indy_gll(n-1)    
enddo
do n = 1,Tdomain%n_elem_z-1
   Tdomain%indz_gll(n) = Tdomain%indz_gll(n)+Tdomain%indz_gll(n-1) 
enddo
!
do n = Tdomain%n_elem_x-1,1,-1
   Tdomain%indx_gll(n) = Tdomain%indx_gll(n-1)        
enddo
do n = Tdomain%n_elem_y-1,1,-1
   Tdomain%indy_gll(n) = Tdomain%indy_gll(n-1)    
enddo
do n = Tdomain%n_elem_z-1,1,-1
   Tdomain%indz_gll(n) = Tdomain%indz_gll(n-1) 
enddo
!
Tdomain%indx_gll(0) = 0        
Tdomain%indy_gll(0) = 0 
Tdomain%indz_gll(0) = 0
!
deallocate(indx_gll_all)
deallocate(indy_gll_all)
deallocate(indz_gll_all)
!
! Done with global numbering common to all proc
!

icount = 0

do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   allocate (Tdomain%specel(n)%Iglobnum(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
   Tdomain%specel(n)%Iglobnum = -1
   do k = 1,ngllz-2
       do j = 1,nglly-2
           do i = 1,ngllx-2
               Tdomain%specel(n)%Iglobnum(i,j,k) = icount
               icount = icount + 1
           enddo
       enddo
   enddo
enddo

do n = 0,Tdomain%n_face-1
   ngllx = Tdomain%sFace(n)%ngll1
   nglly = Tdomain%sFace(n)%ngll2
   allocate (Tdomain%sFace(n)%Iglobnum_Face(1:ngllx-2,1:nglly-2))
   Tdomain%sFace(n)%Iglobnum_Face = -1
   do j = 1,nglly-2
       do i = 1,ngllx-2
           Tdomain%sFace(n)%Iglobnum_Face(i,j) = icount
           icount = icount + 1
       enddo
   enddo
enddo

do n = 0,Tdomain%n_edge-1
    ngllx = Tdomain%sEdge(n)%ngll
    allocate (Tdomain%sEdge(n)%Iglobnum_Edge(1:ngllx-2))
    Tdomain%sEdge(n)%Iglobnum_Edge = -1
    do i = 1,ngllx-2
        Tdomain%sEdge(n)%Iglobnum_Edge(i) = icount
        icount = icount + 1
    enddo
enddo

do n = 0,Tdomain%n_vertex-1
    Tdomain%sVertex(n)%Iglobnum_Vertex = icount
    icount = icount + 1
enddo

Tdomain%n_glob_points = icount 

do n = 0,Tdomain%n_elem - 1
    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz
    ! Taking information from faces
    nf = Tdomain%specel(n)%Near_Faces(0)
    Tdomain%specel(n)%Iglobnum (1:ngllx-2, 1:nglly-2, 0) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngllx-2, 1:nglly-2)
    nf = Tdomain%specel(n)%Near_Faces(5)
    Tdomain%specel(n)%Iglobnum (1:ngllx-2, 1:nglly-2, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngllx-2, 1:nglly-2)
    nf = Tdomain%specel(n)%Near_Faces(1)
    Tdomain%specel(n)%Iglobnum (1:ngllx-2, 0, 1:ngllz-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngllx-2, 1:ngllz-2)
    nf = Tdomain%specel(n)%Near_Faces(3)
    Tdomain%specel(n)%Iglobnum (1:ngllx-2, nglly-1, 1:ngllz-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngllx-2, 1:ngllz-2)
    nf = Tdomain%specel(n)%Near_Faces(2)
    Tdomain%specel(n)%Iglobnum (ngllx-1, 1:nglly-2, 1:ngllz-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:nglly-2, 1:ngllz-2)
    nf = Tdomain%specel(n)%Near_Faces(4)
    Tdomain%specel(n)%Iglobnum (0, 1:nglly-2, 1:ngllz-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:nglly-2, 1:ngllz-2)
    ! Taking information from edges
    ne = Tdomain%specel(n)%Near_Edges(0)
    Tdomain%specel(n)%Iglobnum(1:ngllx-2,0,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    ne = Tdomain%specel(n)%Near_Edges(1)
    Tdomain%specel(n)%Iglobnum(ngllx-1,1:nglly-2,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    ne = Tdomain%specel(n)%Near_Edges(2)
    Tdomain%specel(n)%Iglobnum(1:ngllx-2,nglly-1,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    ne = Tdomain%specel(n)%Near_Edges(3)
    Tdomain%specel(n)%Iglobnum(0,1:nglly-2,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    ne = Tdomain%specel(n)%Near_Edges(4)
    Tdomain%specel(n)%Iglobnum(ngllx-1,0,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    ne = Tdomain%specel(n)%Near_Edges(5)
    Tdomain%specel(n)%Iglobnum(1:ngllx-2,0,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    ne = Tdomain%specel(n)%Near_Edges(6)
    Tdomain%specel(n)%Iglobnum(0,0,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    ne = Tdomain%specel(n)%Near_Edges(7)
    Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    ne = Tdomain%specel(n)%Near_Edges(8)
    Tdomain%specel(n)%Iglobnum(ngllx-1,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    ne = Tdomain%specel(n)%Near_Edges(9)
    Tdomain%specel(n)%Iglobnum(1:ngllx-2,nglly-1,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    ne = Tdomain%specel(n)%Near_Edges(10)
    Tdomain%specel(n)%Iglobnum(0,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    ne = Tdomain%specel(n)%Near_Edges(11)
    Tdomain%specel(n)%Iglobnum(0,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    ! Taking information from vertices
    nv = Tdomain%specel(n)%Near_Vertices(0)
    Tdomain%specel(n)%Iglobnum(0,0,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(1)
    Tdomain%specel(n)%Iglobnum(ngllx-1,0,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(2)
    Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(3)
    Tdomain%specel(n)%Iglobnum(0,nglly-1,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(4)
    Tdomain%specel(n)%Iglobnum(0,0,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(5)
    Tdomain%specel(n)%Iglobnum(ngllx-1,0,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(6)
    Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(7)
    Tdomain%specel(n)%Iglobnum(0,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
enddo

do n = 0,Tdomain%n_face-1
    deallocate (Tdomain%sFace(n)%Iglobnum_Face)
enddo
do n = 0,Tdomain%n_edge-1
    deallocate (Tdomain%sEdge(n)%Iglobnum_Edge)
enddo

return
end
