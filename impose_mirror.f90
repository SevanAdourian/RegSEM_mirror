!=================================================
subroutine impose_mirror (opt, Tdomain, ntime, rg)
!=================================================
!
use sdomains
!
implicit none
!
type (Domain), intent (INOUT), target :: Tdomain
integer, intent (IN) :: ntime, rg
!
integer :: len_DP, n, i, j, x,y,z, ngllx,nglly,ngllz, comp, irec
doubleprecision :: current_t, t_max, wtaper
character*60 :: mirror_file
character :: opt*5
type(mirror), pointer :: mir
doubleprecision, dimension(:), allocatable :: tmp
!
if(opt=='displ')then
   mir=>Tdomain%mirror_displ
elseif(opt=='force')then  
   mir=>Tdomain%mirror_force    
elseif(opt=='excit')then  
   mir=>Tdomain%mirror_excit    
else
   write(*,*)'bad argument opt in impose_mirror.f90'
endif
!
allocate(tmp(0:mir%recl_mirror-1))  
!   
if (ntime==0) then
 if(mir%decim_fact>1)then
    allocate(mir%b((mir%spln_order+1)*mir%decim_fact+1))
    call bmn(mir%b,mir%decim_fact,mir%spln_order) 
    allocate(mir%tmp(0:mir%recl_mirror-1,mir%spln_order+1))
 endif
 inquire(iolength=len_DP)tmp(0:mir%recl_mirror-1)
 if (len_DP/=0) then
   write (mirror_file,"(a,I3.3)") "mirror."//opt//".",rg
   open (mir%lunit,file=mirror_file,access='direct',form="unformatted",recl=len_DP,status='old')
 endif
endif
!
! read spline coefficients
!
if(mir%decim_fact==1)then
 if (mir%recl_mirror/=0)read(mir%lunit,rec=Tdomain%sTimeParam%ntime-ntime)tmp(:)
else
 do j = 1,mir%spln_order+1
    irec = (Tdomain%sTimeParam%ntime-ntime)/mir%decim_fact+1
    print*,rg,mir%decim_fact,Tdomain%sTimeParam%ntime,ntime,irec
    if (mir%recl_mirror/=0)read(mir%lunit,rec=irec+j-1)mir%tmp(:,j)
 enddo
endif

if(mir%decim_fact>1)then
!
! Reconstruct signal
!
tmp(:) = 0.d0
!
do j = 1,mir%spln_order+1
 i = (mir%spln_order+1-j)*mir%decim_fact+mod(Tdomain%sTimeParam%ntime-ntime,mir%decim_fact)+1
 tmp(:) = tmp(:)+mir%tmp(:,j)*mir%b(i)
enddo                    
!
endif
!
i = 0
!
do n = 0,Tdomain%n_elem-1
   if(Tdomain%specel(n)%mirror_position == 1)then
      ngllx = Tdomain%specel(n)%ngllx
      nglly = Tdomain%specel(n)%nglly
      ngllz = Tdomain%specel(n)%ngllz
      do z = 0,ngllz-1
         do y = 0,nglly-1
            do x = 0,ngllx-1
               if (opt=='excit') then
                  do comp = 0,2
                     Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp) = &
                          Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp) &
                          + tmp(i)
                     i = i+1
                  enddo
               endif
               if(Tdomain%specel(n)%win_mirror(x,y,z)/=0)then
                  if (opt=='displ') then
                     do comp = 0,2
                        Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp) = &
                             Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp) & 
                             - tmp(i)*Tdomain%specel(n)%win_mirror(x,y,z)
                        i = i+1
                     enddo
                  elseif (opt=='force') then
                     do comp = 0,2
                        Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp) = &
                             Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp) & 
                             + tmp(i)*Tdomain%specel(n)%win_mirror(x,y,z)
                        i = i+1
                     enddo
                  endif
               endif
            enddo
         enddo
      enddo
   endif
enddo
!
deallocate(tmp)
!
if (i/=mir%recl_mirror)   stop 'A TMP BUFFER HAS A LENGTH DIFFERENT FROM RECL_MIRROR !!!'
!
if (ntime==Tdomain%sTimeParam%ntime-1 .and. mir%recl_mirror/=0)then
 close (mir%lunit)
 if(mir%decim_fact>1)deallocate(mir%tmp)
 if(mir%decim_fact>1)deallocate(mir%b)
endif
!
return
!
!===========================
end subroutine impose_mirror
!===========================
