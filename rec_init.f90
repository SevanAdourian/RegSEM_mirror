!==============================================
subroutine rec_init (Tdomain, ntime, ninit, rg)
!==============================================
!
! Written by Yder MASSON (masson@ipgp.fr) june 20 2012  
! 
! record fields variables for point matching when running the time reversal mirror simulation with attenation 
!
use sdomains
!
implicit none
!
type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: ntime, ninit, rg
!
integer :: len_DP, n,mypos
character*60 :: init_file
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
enddo
!
do n = 0,Tdomain%n_vertex-1
    if (.not.Tdomain%sVertex(n)%PML) then
        len_DP = len_DP + sizeof(Tdomain%svertex(n)%sSimu(0)%Displ)+sizeof(Tdomain%svertex(n)%sSimu(0)%Veloc)
    endif
enddo
!   
   if (len_DP/=0) then
      write (init_file,"(a,I3.3)") "init",rg
      open (45,file=init_file,access='stream',form="unformatted",status='replace')
   endif
!
Tdomain%recl_init = len_DP+sizeof(mypos)
!
endif
!
if(mod(Tdomain%sTimeParam%ntime-ntime-1,ninit)==0)then
!
mypos = (Tdomain%sTimeParam%ntime-ntime-1)/ninit*Tdomain%recl_init+1 
write(45,pos=mypos)Tdomain%recl_init
!
do n = 0,Tdomain%n_elem-1
    if (.not.Tdomain%specel(n)%PML) then
        write (45)Tdomain%specel(n)%sSimu(0)%Displ,-Tdomain%specel(n)%sSimu(0)%Veloc
        if(Tdomain%n_sls>0)then
        write (45)                                 &
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
        write (45)Tdomain%sface(n)%sSimu(0)%Displ,-Tdomain%sface(n)%sSimu(0)%Veloc
    endif
enddo
!
do n = 0,Tdomain%n_edge-1
    if (.not.Tdomain%sEdge(n)%PML) then
        write (45)Tdomain%sedge(n)%sSimu(0)%Displ,-Tdomain%sedge(n)%sSimu(0)%Veloc
    endif
enddo
!
do n = 0,Tdomain%n_vertex-1
    if (.not.Tdomain%sVertex(n)%PML) then
        write (45)Tdomain%svertex(n)%sSimu(0)%Displ,-Tdomain%svertex(n)%sSimu(0)%Veloc
    endif
enddo
!
endif
!
if (ntime==Tdomain%sTimeParam%ntime-1)   close (45)
!
!======================
end subroutine rec_init
!======================
