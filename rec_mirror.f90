!==============================================
subroutine rec_mirror (opt, Tdomain, ntime, rg)
  !==============================================
  !
  use sdomains
  !
  implicit none
  !
  type (Domain), intent (INOUT), target :: Tdomain
  integer, intent (IN) :: ntime, rg
  integer :: len_DP, n, i, j, x,y,z, ngllx,nglly,ngllz, comp
  doubleprecision, dimension(:), allocatable :: tmp
  character*60 :: mirror_file
  character :: opt*5
  character*60 :: meta_file
  type(mirror), pointer :: mir
  !
  if(opt=='displ')then
     mir=>Tdomain%mirror_displ
  elseif(opt=='force')then  
     mir=>Tdomain%mirror_force    
  elseif(opt=='excit')then
     mir=>Tdomain%mirror_excit
  else
     write(*,*)'bad argument opt in rec_mirror.f90', opt
  endif
  !
  ! some init/allocate at first iteration
  !
  allocate(tmp(0:mir%recl_mirror-1))
  !
  if (ntime==0) then
     if (mir%recl_mirror > 0)then
        if(opt=='excit')then
           write(meta_file,"(a,I3.3)") "meta.",rg  
           open(57,file=meta_file,form='formatted',status='replace')
           write(57,*)sngl(Tdomain%sTimeParam%dt)
           write(57,*)mir%decim_fact
           write(57,*)mir%spln_order
           write(57,*)mir%recl_mirror/3
           do i = 1,mir%recl_mirror/3
              write(57,*)Tdomain%index_pos(i)
           enddo
           close(57)
        endif
     endif
     

     if(mir%decim_fact==1)then      
        mir%count_rec = 0 ! init decimate counter  
     else
        allocate(mir%b((mir%spln_order+1)*mir%decim_fact+1))   
        call bmn(mir%b,mir%decim_fact,mir%spln_order)      
        mir%count_rec = 0 ! init decimate counter  
        allocate(mir%tmp(0:mir%recl_mirror-1,mir%spln_order+1))
        mir%tmp(:,:) = 0.d0
     endif
     inquire(iolength=len_DP)tmp(0:mir%recl_mirror-1)
     write (mirror_file,"(a,I3.3)") "mirror."//opt//".",rg  
     if (len_DP/=0) open(mir%lunit,file=mirror_file,access='direct',form="unformatted",recl=len_DP,status='replace')
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
                 if(Tdomain%specel(n)%win_mirror(x,y,z)/=0)then
                    do comp = 0,2
                       if(opt=='excit')then
                          tmp(i) = Tdomain%specel(n)%sSimu(0)%Excitation(x,y,z,comp)
                       else
                          tmp(i) = Tdomain%specel(n)%sSimu(0)%Forces(x,y,z,comp)
                       endif
                       i = i+1     
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  if (i/=mir%recl_mirror)   stop 'A TMP BUFFER HAS A LENGTH DIFFERENT FROM RECL_MIRROR !!!'
  !
  ! filter (convolve with B-spline)
  !
  if(mir%decim_fact>1)then
     do j = 1,mir%spln_order+1
        i = mod(ntime,mir%decim_fact)+(mir%spln_order+1-j)*mir%decim_fact+1
        mir%tmp(:,j) = mir%tmp(:,j)+mir%b(i)*tmp(:)
     enddo
  endif
  !
  ! Decimate/write
  !
  if(mir%decim_fact==1)then
     mir%count_rec = mir%count_rec+1 
     if (mir%recl_mirror/=0) write(mir%lunit,rec=ntime+1)tmp(:)
  else
     if(mod(ntime,mir%decim_fact)==mir%decim_fact-1)then
        mir%count_rec = mir%count_rec+1  
        if (mir%recl_mirror/=0)then
           write(mir%lunit,rec=mir%count_rec)mir%tmp(:,1)
        endif
        do j = 1,mir%spln_order
           mir%tmp(:,j) = mir%tmp(:,j+1)
        enddo
        mir%tmp(:,mir%spln_order+1) = 0
     endif
  endif
  !
  if (ntime==Tdomain%sTimeParam%ntime-1 .and.mir%recl_mirror/=0)then
     if(mir%decim_fact>1)then
        ! dump buffer
        do j = 1,mir%spln_order+1
           mir%count_rec = mir%count_rec+1  
           write (mir%lunit,rec=mir%count_rec)mir%tmp(:,j)
        enddo
        deallocate(mir%tmp) 
        ! post processing : obtain B-spline coefs
        call postproc_bsplines(Tdomain,opt)
        deallocate(mir%b)
     endif
     close (mir%lunit)
     if(opt=='excit')then
        open(57,file=meta_file,form='formatted',status='old',position='append')
        write(57,*)mir%count_rec
        close(57)
     endif
  endif
  !
  deallocate(tmp)
  !
  return
  !
  !========================
end subroutine rec_mirror
!========================
