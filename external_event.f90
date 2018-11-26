subroutine external_event(Tdomain,rg,ntime)

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

  if(ntime==0)then !======== Init ===========!
     ! open file
     if(rg==0)then
        open(49,file='mirror_src.dat',access='stream',form="unformatted",status='old')
        ! read header
        read(49)Tdomain%ssrc_ext%spln_order,Tdomain%ssrc_ext%t0,Tdomain%ssrc_ext%dt,Tdomain%ssrc_ext%n_gll,Tdomain%ssrc_ext%nt
        !write(*,*)'header',Tdomain%ssrc_ext%spln_order,Tdomain%ssrc_ext%t0,Tdomain%ssrc_ext%dt,Tdomain%ssrc_ext%n_gll,Tdomain%ssrc_ext%nt
     endif

     call mpi_bcast(Tdomain%ssrc_ext%spln_order,1,mpi_integer,0,mpi_comm_world,code)
     call mpi_bcast(Tdomain%ssrc_ext%t0        ,1,mpi_real   ,0,mpi_comm_world,code)
     call mpi_bcast(Tdomain%ssrc_ext%dt        ,1,mpi_real   ,0,mpi_comm_world,code)
     call mpi_bcast(Tdomain%ssrc_ext%n_gll     ,1,mpi_integer,0,mpi_comm_world,code)
     call mpi_bcast(Tdomain%ssrc_ext%nt        ,1,mpi_integer,0,mpi_comm_world,code)

     CALL MPI_BARRIER(MPI_COMM_WORLD, code)

     ! allocating and indexing
     !
     Tdomain%ssrc_ext%ngll_mirror = 0
     !
     do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%mirror_position_ext == 1)then
           ngllx = Tdomain%specel(n)%ngllx	
           nglly = Tdomain%specel(n)%nglly
           ngllz = Tdomain%specel(n)%ngllz
           do z = 0,ngllz-1
              do y = 0,nglly-1
                 do x = 0,ngllx-1
                    Tdomain%ssrc_ext%ngll_mirror = Tdomain%ssrc_ext%ngll_mirror+1
                 enddo
              enddo
           enddo
        endif
     enddo
     !
     i = 0
     !      
     allocate(Tdomain%ssrc_ext%glob_num_mirror(Tdomain%ssrc_ext%ngll_mirror))
     !      
     do n = 0,Tdomain%n_elem-1

        if(Tdomain%specel(n)%mirror_position_ext == 1)then

           ngllx = Tdomain%specel(n)%ngllx	
           nglly = Tdomain%specel(n)%nglly
           ngllz = Tdomain%specel(n)%ngllz

           xs = Tdomain%indx_gll(Tdomain%specel(n)%indx)
           ys = Tdomain%indy_gll(Tdomain%specel(n)%indy)
           zs = Tdomain%indz_gll(Tdomain%specel(n)%indz)

           nxtot = Tdomain%ngll_grid_x
           nytot = Tdomain%ngll_grid_y 
           nztot = Tdomain%ngll_grid_z 

           do z = 0,ngllz-1
              do y = 0,nglly-1
                 do x = 0,ngllx-1
                    i = i+1
                    Tdomain%ssrc_ext%glob_num_mirror(i) = (x+xs)+(y+ys)*nxtot+(z+zs)*(nxtot*nytot)
                 enddo
              enddo
           enddo

        endif

     enddo
     !
     if(rg==0)then
        allocate(Tdomain%ssrc_ext%sizes(Tdomain%n_proc))
        allocate(Tdomain%ssrc_ext%displs(Tdomain%n_proc))
     endif
     !
     call mpi_gather(Tdomain%ssrc_ext%ngll_mirror,1,mpi_integer,Tdomain%ssrc_ext%sizes,1,mpi_integer,0,mpi_comm_world,code)  
     !
     if(rg==0)then
        Tdomain%ssrc_ext%ngll_mirror_tot = sum(Tdomain%ssrc_ext%sizes)
        Tdomain%ssrc_ext%displs(:) = Tdomain%ssrc_ext%sizes(:)
        do i = 2,Tdomain%n_proc
           Tdomain%ssrc_ext%displs(i) = Tdomain%ssrc_ext%displs(i)+Tdomain%ssrc_ext%displs(i-1)
        enddo
        do i = Tdomain%n_proc,2,-1
           Tdomain%ssrc_ext%displs(i) = Tdomain%ssrc_ext%displs(i-1)
        enddo
        Tdomain%ssrc_ext%displs(1) = 0
        allocate(Tdomain%ssrc_ext%glob_num_mirror_all(Tdomain%ssrc_ext%ngll_mirror_tot))
     endif

     call mpi_gatherv(Tdomain%ssrc_ext%glob_num_mirror,Tdomain%ssrc_ext%ngll_mirror,mpi_integer,Tdomain%ssrc_ext%glob_num_mirror_all&
          ,Tdomain%ssrc_ext%sizes,Tdomain%ssrc_ext%displs,mpi_integer,0,mpi_comm_world,code)

     !
     !
     !

     if(rg==0)then
        allocate(Tdomain%ssrc_ext%index_num(Tdomain%ssrc_ext%ngll_mirror_tot))
        call indexx_i4b(Tdomain%ssrc_ext%glob_num_mirror_all,Tdomain%ssrc_ext%index_num)
        !  i = index_num(1)    
        !  write(48,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
        !  write(49,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
        !  do n = 2,ngll_mirror_tot
        !    i = index_num(n)
        !    if(glob_num_mirror_all(i)/=glob_num_mirror_all(index_num(n-1)))then
        !    write(48,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
        !    write(49,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
        !    endif   
        !  enddo
        !  write(48,*)
        !  close(49)
        !  deallocate(sizes,displs,x_gll_all,y_gll_all,z_gll_all,glob_num_mirror_all,index_num)
     endif



     allocate(Tdomain%ssrc_ext%b_coefs(Tdomain%n_dim,Tdomain%ssrc_ext%ngll_mirror,Tdomain%ssrc_ext%spln_order+1))
     allocate(Tdomain%ssrc_ext%tmp_x(Tdomain%ssrc_ext%ngll_mirror))
     allocate(Tdomain%ssrc_ext%tmp_y(Tdomain%ssrc_ext%ngll_mirror))
     allocate(Tdomain%ssrc_ext%tmp_z(Tdomain%ssrc_ext%ngll_mirror))
     if(rg==0)allocate(Tdomain%ssrc_ext%tmp_x_all(Tdomain%ssrc_ext%ngll_mirror_tot))
     if(rg==0)allocate(Tdomain%ssrc_ext%tmp_y_all(Tdomain%ssrc_ext%ngll_mirror_tot))
     if(rg==0)allocate(Tdomain%ssrc_ext%tmp_z_all(Tdomain%ssrc_ext%ngll_mirror_tot))
     Tdomain%ssrc_ext%init = .true.


  endif !======== End Init ===========!



  !
  ! read and interpolate
  !
  t = (ntime-1)*Tdomain%sTimeParam%dt
  tnt = (t+Tdomain%ssrc_ext%t0)/Tdomain%ssrc_ext%dt+(Tdomain%ssrc_ext%spln_order+1)*0.5d0
  ind = int(tnt)
  tnt = tnt-dble(ind)
  !
  if(Tdomain%ssrc_ext%init.or.ind/=Tdomain%ssrc_ext%ind_old)then
     !
     j1 = 1
     j2 = Tdomain%ssrc_ext%spln_order+1
     if(.not.Tdomain%ssrc_ext%init)then
        if(ind==Tdomain%ssrc_ext%ind_old+1)then
           Tdomain%ssrc_ext%b_coefs(:,:,1:Tdomain%ssrc_ext%spln_order)=Tdomain%ssrc_ext%b_coefs(:,:,2:Tdomain%ssrc_ext%spln_order+1)
           j1 = Tdomain%ssrc_ext%spln_order+1
        endif
        if(ind==Tdomain%ssrc_ext%ind_old-1)then
           Tdomain%ssrc_ext%b_coefs(:,:,2:Tdomain%ssrc_ext%spln_order+1)=Tdomain%ssrc_ext%b_coefs(:,:,1:Tdomain%ssrc_ext%spln_order)
           j2 = 1
        endif
     endif
     !
     do j = j1,j2

        if(ind+j<1.or.ind+j>Tdomain%ssrc_ext%nt)then
           Tdomain%ssrc_ext%b_coefs(1,:,j) = 0.d0
           Tdomain%ssrc_ext%b_coefs(2,:,j) = 0.d0
           Tdomain%ssrc_ext%b_coefs(3,:,j) = 0.d0
        else
           ! read
           if(rg==0)then
              tRot = transpose(Tdomain%rot)
              i = Tdomain%ssrc_ext%index_num(1)    
              read(49,pos=INT8(21)+(INT8(ind)+INT8(j)-INT8(1))*INT8(Tdomain%ssrc_ext%n_gll)*INT8(3)*INT8(4))ux,uy,uz
              Tdomain%ssrc_ext%tmp_x_all(i) = tRot(0,0)*ux + tRot(0,1)*uy + tRot(0,2)*uz
              Tdomain%ssrc_ext%tmp_y_all(i) = tRot(1,0)*ux + tRot(1,1)*uy + tRot(1,2)*uz
              Tdomain%ssrc_ext%tmp_z_all(i) = tRot(2,0)*ux + tRot(2,1)*uy + tRot(2,2)*uz
              k = 1
              do n = 2,Tdomain%ssrc_ext%ngll_mirror_tot
                 i = Tdomain%ssrc_ext%index_num(n)
                 if(Tdomain%ssrc_ext%glob_num_mirror_all(i)/=Tdomain%ssrc_ext%glob_num_mirror_all(Tdomain%ssrc_ext%index_num(n-1)))then
                    read(49)ux,uy,uz
                    k = k+1
                 endif
                 Tdomain%ssrc_ext%tmp_x_all(i) = tRot(0,0)*ux + tRot(0,1)*uy + tRot(0,2)*uz
                 Tdomain%ssrc_ext%tmp_y_all(i) = tRot(1,0)*ux + tRot(1,1)*uy + tRot(1,2)*uz
                 Tdomain%ssrc_ext%tmp_z_all(i) = tRot(2,0)*ux + tRot(2,1)*uy + tRot(2,2)*uz
              enddo

              if(k/=Tdomain%ssrc_ext%n_gll)then
                 write(*,*)'length doesn"t match',k,Tdomain%ssrc_ext%n_gll
              else
                 !write(*,*)' it seems ok we have',maxval(Tdomain%ssrc_ext%tmp_x_all(:)),maxval(Tdomain%ssrc_ext%tmp_x_all(:)),ind+j-1
              endif
           endif
           ! distribute   
           call mpi_scatterv(Tdomain%ssrc_ext%tmp_x_all,Tdomain%ssrc_ext%sizes,Tdomain%ssrc_ext%displs,mpi_real,&
                Tdomain%ssrc_ext%tmp_x,Tdomain%ssrc_ext%ngll_mirror,mpi_real,0,mpi_comm_world,code)
           call mpi_scatterv(Tdomain%ssrc_ext%tmp_y_all,Tdomain%ssrc_ext%sizes,Tdomain%ssrc_ext%displs,mpi_real,&
                Tdomain%ssrc_ext%tmp_y,Tdomain%ssrc_ext%ngll_mirror,mpi_real,0,mpi_comm_world,code)
           call mpi_scatterv(Tdomain%ssrc_ext%tmp_z_all,Tdomain%ssrc_ext%sizes,Tdomain%ssrc_ext%displs,mpi_real,&
                Tdomain%ssrc_ext%tmp_z,Tdomain%ssrc_ext%ngll_mirror,mpi_real,0,mpi_comm_world,code)

           Tdomain%ssrc_ext%b_coefs(1,:,j) = Tdomain%ssrc_ext%tmp_x(:)
           Tdomain%ssrc_ext%b_coefs(2,:,j) = Tdomain%ssrc_ext%tmp_y(:)
           Tdomain%ssrc_ext%b_coefs(3,:,j) = Tdomain%ssrc_ext%tmp_z(:)

        endif


     enddo

     Tdomain%ssrc_ext%ind_old = ind

  endif
  !
  ! end read 
  !                  
  ! interpolate

  Tdomain%ssrc_ext%tmp_x = 0
  Tdomain%ssrc_ext%tmp_y = 0
  Tdomain%ssrc_ext%tmp_z = 0

  do j = 0,Tdomain%ssrc_ext%spln_order
     Tdomain%ssrc_ext%tmp_x(:) = Tdomain%ssrc_ext%tmp_x(:)+bspln(0,Tdomain%ssrc_ext%spln_order,dble(tnt+j))*Tdomain%ssrc_ext%b_coefs(1,:,Tdomain%ssrc_ext%spln_order+1-j)
     Tdomain%ssrc_ext%tmp_y(:) = Tdomain%ssrc_ext%tmp_y(:)+bspln(0,Tdomain%ssrc_ext%spln_order,dble(tnt+j))*Tdomain%ssrc_ext%b_coefs(2,:,Tdomain%ssrc_ext%spln_order+1-j)
     Tdomain%ssrc_ext%tmp_z(:) = Tdomain%ssrc_ext%tmp_z(:)+bspln(0,Tdomain%ssrc_ext%spln_order,dble(tnt+j))*Tdomain%ssrc_ext%b_coefs(3,:,Tdomain%ssrc_ext%spln_order+1-j)
  enddo

  ! impose source

  ! Internal Forces

  i_simu = 0
  i = 0
  do n = 0,Tdomain%n_elem-1
     mat = Tdomain%specel(n)%mat_index
     if(Tdomain%specel(n)%mirror_position_ext == 1)then
        ngllx = Tdomain%specel(n)%ngllx	
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        allocate(sav_forces(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
        sav_forces(:,:,:,:) = Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,:)      

        ! 1 fill force with displ
        do z = 0,ngllz-1
           do y = 0,nglly-1
              do x = 0,ngllx-1
                 i = i+1
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,0) = Tdomain%ssrc_ext%tmp_x(i)
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,1) = Tdomain%ssrc_ext%tmp_y(i)
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,2) = Tdomain%ssrc_ext%tmp_z(i)
              enddo
           enddo
        enddo
        i = i-ngllx*nglly*ngllz
        if (.not. Tdomain%specel(n)%PML) then
           call internal_forces (Tdomain%specel(n), i_simu, Tdomain%sSubDomain(mat)%hprimex, &
                Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls, Tdomain%aniso, &
                .false., Tdomain%t_reversal_mirror,.true.)
        else
           write(*,*)' we have a probleme with mirror sources in the pmls...'
        endif

        ! 2 add force to force
        sav_forces(:,:,:,0) = sav_forces(:,:,:,0) + Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,0)*Tdomain%specel(n)%win_mirror_ext(:,:,:)
        sav_forces(:,:,:,1) = sav_forces(:,:,:,1) + Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,1)*Tdomain%specel(n)%win_mirror_ext(:,:,:)
        sav_forces(:,:,:,2) = sav_forces(:,:,:,2) + Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,2)*Tdomain%specel(n)%win_mirror_ext(:,:,:)
        ! 3 fill force with windowed displ
        do z = 0,ngllz-1
           do y = 0,nglly-1
              do x = 0,ngllx-1
                 i = i+1
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,0) = Tdomain%ssrc_ext%tmp_x(i)*Tdomain%specel(n)%win_mirror_ext(x,y,z)
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,1) = Tdomain%ssrc_ext%tmp_y(i)*Tdomain%specel(n)%win_mirror_ext(x,y,z)
                 Tdomain%specel(n)%sSimu(i_simu)%Forces(x,y,z,2) = Tdomain%ssrc_ext%tmp_z(i)*Tdomain%specel(n)%win_mirror_ext(x,y,z)
              enddo
           enddo
        enddo
        if (.not. Tdomain%specel(n)%PML) then
           call internal_forces (Tdomain%specel(n), i_simu, Tdomain%sSubDomain(mat)%hprimex, &
                Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls, Tdomain%aniso, &
                .false., Tdomain%t_reversal_mirror,.true.)
        endif
        ! 4 add force to force     
        sav_forces(:,:,:,:) = sav_forces(:,:,:,:) - Tdomain%specel(n)%sSimu(i_simu)%Forces(:,:,:,:)
        
        deallocate(sav_forces)
     endif
  enddo
  if(i/=Tdomain%ssrc_ext%ngll_mirror)write(*,*)'we have a problem with in in external event...'
  Tdomain%ssrc_ext%init = .false.

  if(ntime==Tdomain%sTimeParam%ntime-1)then !======== Final ===========!
     ! close file
     close(49)
     ! deallocate
     deallocate(Tdomain%ssrc_ext%glob_num_mirror)
  endif !======== End Final ===========!

  return
end subroutine external_event
