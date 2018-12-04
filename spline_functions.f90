
!========================================
recursive function bspln(i,k,x) result(b)
!========================================

implicit none

double precision :: x,b
integer :: k,i

b=0.d0
if(k+1>0)then
   if(k+1==1)then
      if(x>=i.and.x<i+1)b=1.d0   
   else
      b=(x-i)*bspln(i,k-1,x)/(k+1-1)+(i+k+1-x)*bspln(i+1,k-1,x)/(k+1-1)
      if(k==0)b=0
   endif
endif

!=================
end function bspln
!=================

!====================
subroutine bmn(b,m,n)
!====================
  !
  implicit none
  !
  integer :: m,n,i
  double precision :: b((n+1)*m+1),bspln
  !
  do i = 1,(n+1)*m+1
     b(i) = bspln(0,n,dble(i-1)/dble((n+1)*m)*dble(n+1))
  enddo
  !
return
!=================
end subroutine bmn
!=================

!==========================================
subroutine bchfac ( w, nbands, nrow, diag )
!==========================================

!*****************************************************************************80
!
!! BCHFAC constructs a Cholesky factorization of a matrix.
!
!  Discussion:
!
!    The factorization has the form
!
!      C = L * D * L'
!  
!    with L unit lower triangular and D diagonal, for a given matrix C of 
!    order NROW, where C is symmetric positive semidefinite and banded, 
!    having NBANDS diagonals at and below the main diagonal.
! 
!    Gauss elimination is used, adapted to the symmetry and bandedness of C.
! 
!    Near-zero pivots are handled in a special way.  The diagonal 
!    element C(N,N) = W(1,N) is saved initially in DIAG(N), all N. 
! 
!    At the N-th elimination step, the current pivot element, W(1,N), 
!    is compared with its original value, DIAG(N).  If, as the result 
!    of prior elimination steps, this element has been reduced by about 
!    a word length, that is, if W(1,N) + DIAG(N) <= DIAG(N), then the pivot 
!    is declared to be zero, and the entire N-th row is declared to
!    be linearly dependent on the preceding rows.  This has the effect 
!    of producing X(N) = 0 when solving C * X = B for X, regardless of B.
! 
!    Justification for this is as follows.  In contemplated applications 
!    of this program, the given equations are the normal equations for 
!    some least-squares approximation problem, DIAG(N) = C(N,N) gives 
!    the norm-square of the N-th basis function, and, at this point, 
!    W(1,N) contains the norm-square of the error in the least-squares 
!    approximation to the N-th basis function by linear combinations 
!    of the first N-1.  
!
!    Having W(1,N)+DIAG(N) <= DIAG(N) signifies that the N-th function 
!    is linearly dependent to machine accuracy on the first N-1 
!    functions, therefore can safely be left out from the basis of 
!    approximating functions.
!
!    The solution of a linear system C * X = B is effected by the 
!    succession of the following two calls:
! 
!      call bchfac ( w, nbands, nrow, diag )
!
!      call bchslv ( w, nbands, nrow, b, x )
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) W(NBANDS,NROW).
!    On input, W contains the NBANDS diagonals in its rows, 
!    with the main diagonal in row 1.  Precisely, W(I,J) 
!    contains C(I+J-1,J), I=1,...,NBANDS, J=1,...,NROW.
!    For example, the interesting entries of a seven diagonal
!    symmetric matrix C of order 9 would be stored in W as
! 
!      11 22 33 44 55 66 77 88 99
!      21 32 43 54 65 76 87 98  *
!      31 42 53 64 75 86 97  *  *
!      41 52 63 74 85 96  *  *  *
!
!    Entries of the array not associated with an
!    entry of C are never referenced.
!    On output, W contains the Cholesky factorization 
!    C = L*D*L', with W(1,I) containing 1/D(I,I) and W(I,J) 
!    containing L(I-1+J,J), I=2,...,NBANDS.
!
!    Input, integer ( kind = 4 ) NBANDS, indicates the bandwidth of the
!    matrix C, that is, C(I,J) = 0 for NBANDS < abs(I-J).
! 
!    Input, integer ( kind = 4 ) NROW, is the order of the matrix C.
! 
!    Work array, real ( kind = 8 ) DIAG(NROW).
!
  implicit none

  integer ( kind = 4 ) nbands
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) diag(nrow)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) n
  real ( kind = 8 ) ratio
  real ( kind = 8 ) w(nbands,nrow)

  if ( nrow <= 1 ) then
    if ( 0.0D+00 < w(1,1) ) then
      w(1,1) = 1.0D+00 / w(1,1)
    end if
    return
  end if
!
!  Store the diagonal.
!
  diag(1:nrow) = w(1,1:nrow)
!
!  Factorization.
!
  do n = 1, nrow
 
    if ( w(1,n) + diag(n) <= diag(n) ) then
      w(1:nbands,n) = 0.0D+00
    else
 
      w(1,n) = 1.0D+00 / w(1,n)
 
      imax = min ( nbands - 1, nrow - n )
 
      jmax = imax
 
      do i = 1, imax
 
        ratio = w(i+1,n) * w(1,n)
 
        do j = 1, jmax
          w(j,n+i) = w(j,n+i) - w(j+i,n) * ratio
        end do
 
        jmax = jmax-1
        w(i+1,n) = ratio
 
      end do
 
    end if
 
  end do
 
  return

!==
end
!==

!================================================
subroutine bchslv (Tdomain, w, nbands, nrow, opt)
!================================================

!*****************************************************************************
!
!! BCHSLV solves a banded symmetric positive definite system.
!
!  Discussion:
!
!    The system is of the form:
!
!      C * X = B 
!  
!    and the Cholesky factorization of C has been constructed 
!    by BCHFAC.
! 
!    With the factorization 
!
!      C = L * D * L'
!
!    available, where L is unit lower triangular and D is diagonal, 
!    the triangular system 
!
!      L * Y = B 
!
!    is solved for Y (forward substitution), Y is stored in B, the 
!    vector D**(-1)*Y is computed and stored in B, then the 
!    triangular system L'*X = D**(-1)*Y is solved for X 
!    (back substitution).
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W(NBANDS,NROW), the Cholesky factorization for C, 
!    as computed by BCHFAC.
! 
!    Input, integer ( kind = 4 ) NBANDS, the bandwidth of C.
!
!    Input, integer ( kind = 4 ) NROW, the order of the matrix C.
! 
!    Input/output, real ( kind = 8 ) B(NROW).
!    On input, the right hand side.
!    On output, the solution.
!

  use sdomains
  
  implicit none

  type (Domain), intent (INOUT), target :: Tdomain        
  type(mirror), pointer :: mir
  integer :: j
  integer :: n
  integer :: nrow,nbands,lunit
  double precision :: w(nbands,nrow)
  character :: opt*5      
! 
!nbands = mir%spln_order+1, nrow =  mir%count_rec
!
  if(opt=='displ')then
   mir => Tdomain%mirror_displ
  elseif(opt=='force')then
   mir => Tdomain%mirror_force
  endif
!
  lunit = mir%lunit
!
  allocate(mir%tmp(0:mir%recl_mirror-1,0:mir%spln_order))
!
  if ( nrow <= 1 ) then
    write(*,*)'warning few time steps after decimation'
    return
  end if
!
!  Forward substitution. 
!  Solve L*Y = B.
!
  do n = 1, nrow-nbands+1
     if(n==1)then
        do j = 0, nbands - 1
           read(lunit,rec=j+n)mir%tmp(:,j)
        enddo
     else
        do j = 0,nbands-2
           mir%tmp(:,j) = mir%tmp(:,j+1)
        enddo
         read(lunit,rec=n+nbands-1)mir%tmp(:,nbands-1)
     endif
     do j = 1, nbands - 1
        mir%tmp(:,j) = mir%tmp(:,j) - w(j+1,n) * mir%tmp(:,0)
     end do
     
     write(lunit,rec=n)mir%tmp(:,0)
    
 end do

  do n = nrow-nbands+2, nrow

     do j = 0,nbands-2
        mir%tmp(:,j) = mir%tmp(:,j+1)
     enddo

     do j = 1, nrow - n
        mir%tmp(:,j) = mir%tmp(:,j) - w(j+1,n) * mir%tmp(:,0)
     end do

      write(mir%lunit,rec=n)mir%tmp(:,0)
     
  end do
!
!  Back substitution. 
!  Solve L'*X = D**(-1)*Y.
!
  do n = nrow, nrow-nbands+2, -1

     do j = nrow-n,1,-1
        mir%tmp(:,j) = mir%tmp(:,j-1)
     enddo
     read(lunit,rec=n)mir%tmp(:,0)
     
    mir%tmp(:,0) = mir%tmp(:,0) * w(1,n)

    do j = 1, nrow - n
      mir%tmp(:,0) = mir%tmp(:,0) - w(j+1,n) * mir%tmp(:,j)
    end do

     write(lunit,rec=n)mir%tmp(:,0)

  end do

  do n = nrow-nbands+1, 1, -1

     do j = nbands-1,1,-1
        mir%tmp(:,j) = mir%tmp(:,j-1)
     enddo
     read(lunit,rec=n)mir%tmp(:,0)
     
    mir%tmp(:,0) = mir%tmp(:,0) * w(1,n)

    do j = 1, nbands - 1
      mir%tmp(:,0) = mir%tmp(:,0) - w(j+1,n) * mir%tmp(:,j)
    end do

     write(lunit,rec=n)mir%tmp(:,0)

  end do

  deallocate(mir%tmp)
  return

!====================
end subroutine bchslv
!====================

!========================================
subroutine postproc_bsplines(Tdomain,opt)
!========================================
!
use sdomains
!
implicit none
!
type (Domain), intent (INOUT), target :: Tdomain
type (mirror), pointer :: mir
integer, pointer :: m,n,np,lunit,nbs
double precision, allocatable :: diag(:,:),ddiag(:)
character :: opt*5
!
if(opt=='displ')then
mir => Tdomain%mirror_displ
elseif(opt=='force')then
mir => Tdomain%mirror_force
else
write(*,*)'bad value for argument "opt" in subroutine postproc_bsplines'
endif   
!
allocate(diag(mir%spln_order+1,mir%count_rec))
allocate(ddiag(mir%count_rec))
call fill_spmat_diag(mir%decim_fact,mir%spln_order,mir%count_rec,Tdomain%sTimeParam%ntime,diag)
call bchfac ( diag, mir%spln_order+1, mir%count_rec, ddiag )
call bchslv (Tdomain, diag, mir%spln_order+1, mir%count_rec, opt)
deallocate(diag)
deallocate(ddiag)
!
return
!  
!===============================
end subroutine postproc_bsplines
!===============================

!===========================================
subroutine fill_spmat_diag(m,n,nsp,nt,spmat)
!===========================================
!
implicit none
!
integer :: m,n,nsp,nt
double precision :: spmat(n+1,nsp),b((n+1)*m+1)
!
integer :: i,j,i1,i2,i3,i4
!
call bmn(b,m,n)
!
do j = 1,nsp
   do i = 1,n+1
!
      i1 = max(1+(i-1)*m,(n+1-j)*m+1)
      i2 = min(1+(n+1)*m,(n-j+1)*m+nt)
      i3 = max(1,(n+2-j-i)*m+1)
      i4 = min(1+(n+2-i)*m,(n+2-j-i)*m+nt)
!
      if(j+i-1<=nsp)then
         spmat(i,j) = dot_product(b(i1:i2),b(i3:i4))
      else
         spmat(i,j) = 0.d0
      endif
!
   enddo
enddo
!
return
!=============================
end subroutine fill_spmat_diag
!=============================
