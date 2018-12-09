subroutine get_mirror_exc(Tdomain,rg,ntime)

  use sdomains
  use sort

  implicit none

  include 'mpif.h'

  integer :: rg,ntime
  type (domain), intent (INOUT) :: Tdomain

  integer :: code,x,y,z
  integer :: ngllx,nglly,ngllz,xs,ys,zs,nxtot,nytot,nztot,i_simu,mat
  real :: ux,uy,uz,t,tnt
  integer :: ind,i,j,k,n,j1,j2
  doubleprecision, dimension(0:2,0:2) :: tRot
  double precision :: bspln
  double precision, allocatable :: sav_forces(:,:,:,:)


  ! impose source

  ! Internal Forces

  i_simu = 0
  i = 0
  do n = 0,Tdomain%n_elem-1
     mat = Tdomain%specel(n)%mat_index
     if(Tdomain%specel(n)%mirror_position == 1)then
        ngllx = Tdomain%specel(n)%ngllx	
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        allocate(sav_forces(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
        sav_forces(:,:,:,:) = Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,:)
        ! Tdomain%specel(n)%sSimu(i_simu)%Excitation(:,:,:,:) = Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,:)
        ! 1 fill force with displ
        i = i-ngllx*nglly*ngllz
        if (.not. Tdomain%specel(n)%PML) then
           call internal_forces (Tdomain%specel(n), i_simu, Tdomain%sSubDomain(mat)%hprimex, &
                Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls, Tdomain%aniso, &
                .false., Tdomain%t_reversal_mirror,.false.)
        else
           write(*,*)' we have a probleme with mirror sources in the pmls...'
        endif

        ! 2 add force to force
        Tdomain%specel(n)%sSimu(i_simu)%Excitation(:,:,:,0) = - Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,0)*Tdomain%specel(n)%win_mirror(:,:,:)
        Tdomain%specel(n)%sSimu(i_simu)%Excitation(:,:,:,1) = - Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,1)*Tdomain%specel(n)%win_mirror(:,:,:)
        Tdomain%specel(n)%sSimu(i_simu)%Excitation(:,:,:,2) = - Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,2)*Tdomain%specel(n)%win_mirror(:,:,:)

        ! 3 fill force with windowed displ
        do z = 0,ngllz-1
           do y = 0,nglly-1
              do x = 0,ngllx-1
                 i = i+1
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,0) = sav_forces(x,y,z,0)*Tdomain%specel(n)%win_mirror(x,y,z)
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,1) = sav_forces(x,y,z,1)*Tdomain%specel(n)%win_mirror(x,y,z)
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,2) = sav_forces(x,y,z,2)*Tdomain%specel(n)%win_mirror(x,y,z)
              enddo
           enddo
        enddo
        if (.not. Tdomain%specel(n)%PML) then
           call internal_forces (Tdomain%specel(n), i_simu, Tdomain%sSubDomain(mat)%hprimex, &
                Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls, Tdomain%aniso, &
                .false., Tdomain%t_reversal_mirror,.false.)
        endif
        ! 4 add force to force     
        Tdomain%specel(n)%sSimu(i_simu)%Excitation(:,:,:,:) = Tdomain%specel(n)%sSimu(i_simu)%Excitation(:,:,:,:) + Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,:)

        ! Back to the original forces 
        Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,:) = sav_forces(:,:,:,:)
        deallocate(sav_forces)
     endif
  enddo
  ! if(i/=Tdomain%ssrc_ext%ngll_mirror)write(*,*)'we have a problem with in in external event...'
  ! if (Tdomain%t_reversal_mirror==1)   call rec_mirror ('exc', Tdomain, ntime, rg)
  ! Tdomain%ssrc_ext%init = .false.

  return
end subroutine get_mirror_exc
