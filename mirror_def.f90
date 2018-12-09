subroutine mirror_definition (Tdomain,rg)


  use sdomains
  use read_mirror
  use sort

  implicit none

  include 'mpif.h'

  type (Domain), intent (INOUT) :: Tdomain

  integer :: n, i, ii, i1, i2, j1, j2, k1, k2, x, y, z, x1, x2, y1, y2, z1, z2, ngll_mirror, ngll_inner, ngll_if
  integer :: indx, indy, indz, ngllx, nglly, ngllz, rg, code
  integer :: xs, ys, zs, nxtot, nytot, nztot, iglob, ngll_mirror_tot, count_gll
  integer, dimension(:), allocatable :: glob_num_mirror, glob_num_mirror_all, sizes, displs
  integer, dimension(:), allocatable :: index_num, index_pos_all, index_pos
  double precision, dimension(:), allocatable :: x_gll, y_gll, z_gll, x_gll_all, y_gll_all, z_gll_all
  double precision :: xa, ya, za, winx, winy, winz, sum1, sum2
  character*60 :: temp_file
  !
  ! some parameters for recording/imposing the mirror...
  !
  Tdomain%mirror_displ%lunit = 44 ! logical file units
  Tdomain%mirror_force%lunit = 43 !
  Tdomain%mirror_excit%lunit = 42 ! logical file units
  Tdomain%mirror_displ%lunit2 = 46 ! logical file units
  Tdomain%mirror_force%lunit2 = 45 !
  Tdomain%mirror_excit%lunit2 = 47 ! logical file units
  Tdomain%mirror_displ%spln_order = 1 ! B-Spline order
  Tdomain%mirror_force%spln_order = 1 !
  Tdomain%mirror_excit%spln_order = 1 !
  !
  ! Mirror geometry...
  !
  i1 = Tdomain%i1_mirror
  i2 = Tdomain%n_elem_x-Tdomain%i2_mirror-1
  j1 = Tdomain%j1_mirror
  j2 = Tdomain%n_elem_y-Tdomain%j2_mirror-1
  k1 = Tdomain%k1_mirror
  k2 = Tdomain%n_elem_z-Tdomain%k2_mirror-1

  do n = 0,Tdomain%n_elem-1

     indx = Tdomain%specel(n)%indx
     indy = Tdomain%specel(n)%indy
     indz = Tdomain%specel(n)%indz

     ngllx = Tdomain%specel(n)%ngllx
     nglly = Tdomain%specel(n)%nglly
     ngllz = Tdomain%specel(n)%ngllz

     allocate (Tdomain%specel(n)%win_mirror(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
     allocate (Tdomain%specel(n)%win_mirror_ext(0:ngllx-1, 0:nglly-1, 0:ngllz-1))

     Tdomain%specel(n)%win_mirror = 0
     Tdomain%specel(n)%win_mirror_ext = 0

     ! inner mirror for reconstructing time reversed wavefield

     if(indx>=i1+1 .and. indx<=i2-1 .and. indy>=j1+1 .and. indy<=j2-1 .and. indz>=k1+1 .and. indz<=k2-1)then

        x1 = 0
        x2 = ngllx-1
        y1 = 0
        y2 = nglly-1 
        z1 = 0
        z2 = ngllz-1

        !if(indx==i1+1)x1=x2
        !if(indx==i2-1)x2=x1
        !if(indy==j1+1)y1=y2
        !if(indy==j2-1)y2=y1
        !if(indz==k1+1)z1=z2
        !if(indz==k2-1)z2=z1


        do z = z1,z2
           winz = 1
           if(indz==k1+1)winz = 1-(0.5*cos(4.*atan(1.)*real(z-z1)/real(z2-z1))+0.5)
           if(indz==k2-1)winz =   (0.5*cos(4.*atan(1.)*real(z-z1)/real(z2-z1))+0.5)
           do y = y1,y2
              winy = 1
              if(indy==j1+1)winy = 1-(0.5*cos(4.*atan(1.)*real(y-y1)/real(y2-y1))+0.5)
              if(indy==j2-1)winy =   (0.5*cos(4.*atan(1.)*real(y-y1)/real(y2-y1))+0.5)
              do x = x1,x2
                 winx = 1
                 if(indx==i1+1)winx = 1-(0.5*cos(4.*atan(1.)*real(x-x1)/real(x2-x1))+0.5)
                 if(indx==i2-1)winx =   (0.5*cos(4.*atan(1.)*real(x-x1)/real(x2-x1))+0.5)
                 Tdomain%specel(n)%win_mirror(x,y,z)=winx*winy*winz
              enddo
           enddo
        enddo

     endif

     ! outer mirror for reconstructing global wavefield associated with a teleseismic event

     if(indx>=i1 .and. indx<=i2 .and. indy>=j1 .and. indy<=j2 .and. indz>=k1 .and. indz<=k2)then

        x1 = 0
        x2 = ngllx-1
        y1 = 0
        y2 = nglly-1 
        z1 = 0
        z2 = ngllz-1

        !if(indx==i1)x1=x2
        !if(indx==i2)x2=x1
        !if(indy==j1)y1=y2
        !if(indy==j2)y2=y1
        !if(indz==k1)z1=z2
        !if(indz==k2)z2=z1


        do z = z1,z2
           winz = 1
           if(indz==k1)winz = 1-(0.5*cos(4.*atan(1.)*real(z-z1)/real(z2-z1))+0.5)
           if(indz==k2)winz =   (0.5*cos(4.*atan(1.)*real(z-z1)/real(z2-z1))+0.5)
           do y = y1,y2
              winy = 1
              if(indy==j1)winy = 1-(0.5*cos(4.*atan(1.)*real(y-y1)/real(y2-y1))+0.5)
              if(indy==j2)winy =   (0.5*cos(4.*atan(1.)*real(y-y1)/real(y2-y1))+0.5)
              do x = x1,x2
                 winx = 1
                 if(indx==i1)winx = 1-(0.5*cos(4.*atan(1.)*real(x-x1)/real(x2-x1))+0.5)
                 if(indx==i2)winx =   (0.5*cos(4.*atan(1.)*real(x-x1)/real(x2-x1))+0.5)
                 Tdomain%specel(n)%win_mirror_ext(x,y,z)=winx*winy*winz
              enddo
           enddo
        enddo

     endif

  enddo

  do n = 0,Tdomain%n_elem-1

     ngllx = Tdomain%specel(n)%ngllx
     nglly = Tdomain%specel(n)%nglly
     ngllz = Tdomain%specel(n)%ngllz

     ! For a mirror of arbitrary shape, feed the code with subroutine get_mirror and get_mirror_ext(these must be placed in the /SOLVE/src/Modules/read_mirror.f90)

     ! do z = 0,ngllz-1
     !    do y = 0,nglly-1
     !       do x = 0,ngllx-1
     !           i = Tdomain%specel(n)%Iglobnum(x,y,z)
     !          call get_mirror(Tdomain%Globcoord(0,i),Tdomain%Globcoord(1,i),Tdomain%Globcoord(2,i),&
     !                          Tdomain%specel(n)%win_mirror(x,y,z))
     !          call get_mirror_ext(Tdomain%Globcoord(0,i),Tdomain%Globcoord(1,i),Tdomain%Globcoord(2,i),&
     !                              Tdomain%specel(n)%win_mirror(x,y,z))
     !       enddo
     !    enddo
     ! enddo

  enddo

  call mpi_barrier(mpi_comm_world, code)

  do n = 0,Tdomain%n_elem-1

     ngllx = Tdomain%specel(n)%ngllx	
     nglly = Tdomain%specel(n)%nglly
     ngllz = Tdomain%specel(n)%ngllz

     sum1 = 0
     sum2 = 0

     do z = 0,ngllz-1
        do y = 0,nglly-1
           do x = 0,ngllx-1
              if(Tdomain%specel(n)%win_mirror(x,y,z)/=0)then
                 sum1 = sum1+1
              endif
              if(Tdomain%specel(n)%win_mirror_ext(x,y,z)/=0)then
                 sum2 = sum2+1
              endif
           enddo
        enddo
     enddo

     if(sum1==0)then
        Tdomain%specel(n)%mirror_position = 0
        deallocate (Tdomain%specel(n)%win_mirror)
     elseif(sum1==ngllx*nglly*ngllz)then
        Tdomain%specel(n)%mirror_position = 2
        deallocate (Tdomain%specel(n)%win_mirror)
     else
        Tdomain%specel(n)%mirror_position = 1
     endif


     if(sum2==0)then
        Tdomain%specel(n)%mirror_position_ext = 0
        deallocate (Tdomain%specel(n)%win_mirror_ext)
     elseif(sum2==ngllx*nglly*ngllz)then
        Tdomain%specel(n)%mirror_position_ext = 2
        deallocate (Tdomain%specel(n)%win_mirror_ext)
     else
        Tdomain%specel(n)%mirror_position_ext = 1
     endif

  enddo

  ngll_inner = 0
  ngll_mirror = 0
  ngll_if = 0

  call mpi_barrier(mpi_comm_world, code)

  do n = 0,Tdomain%n_elem-1

     if(Tdomain%specel(n)%mirror_position == 1)then

        ngllx = Tdomain%specel(n)%ngllx	
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        do z = 0,ngllz-1
           do y = 0,nglly-1
              do x = 0,ngllx-1
                 ngll_inner = ngll_inner + 1
                 if(Tdomain%specel(n)%win_mirror(x,y,z)/=0) ngll_if = ngll_if + 1
              enddo
           enddo
        enddo

     endif
     
     if(Tdomain%specel(n)%mirror_position_ext == 1)then

        ngllx = Tdomain%specel(n)%ngllx	
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        do z = 0,ngllz-1
           do y = 0,nglly-1
              do x = 0,ngllx-1
                 ngll_mirror = ngll_mirror+1
              enddo
           enddo
        enddo

     endif

  enddo
  !
  i = 3*ngll_if
  Tdomain%mirror_displ%recl_mirror = i
  Tdomain%mirror_force%recl_mirror = i
  Tdomain%mirror_excit%recl_mirror = 3*ngll_inner ! We need all GLL points
  !
  call mpi_barrier(mpi_comm_world, code)

! ==================================================================================================

  if(Tdomain%t_reversal_mirror==3)then ! Output outer mirror's GLLs for globalwave propagation solver 
     !
     allocate(glob_num_mirror(ngll_mirror))
     allocate(x_gll(ngll_mirror))
     allocate(y_gll(ngll_mirror))
     allocate(z_gll(ngll_mirror))
     !
     i = 0
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
                    iglob = Tdomain%specel(n)%Iglobnum(x,y,z)

                    glob_num_mirror(i) = (x+xs)+(y+ys)*nxtot+(z+zs)*(nxtot*nytot)
                    xa = Tdomain%Globcoord(0,iglob)
                    ya = Tdomain%Globcoord(1,iglob)
                    za = Tdomain%Globcoord(2,iglob)
                    x_gll(i) = Tdomain%rot(0,0)*xa + Tdomain%rot(0,1)*ya + Tdomain%rot(0,2)*za       
                    y_gll(i) = Tdomain%rot(1,0)*xa + Tdomain%rot(1,1)*ya + Tdomain%rot(1,2)*za
                    z_gll(i) = Tdomain%rot(2,0)*xa + Tdomain%rot(2,1)*ya + Tdomain%rot(2,2)*za

                 enddo
              enddo
           enddo

        endif

     enddo

     if(rg==0)then
        allocate(sizes(Tdomain%n_proc))
        allocate(displs(Tdomain%n_proc))
     endif

     call mpi_gather(ngll_mirror,1,mpi_integer,sizes,1,mpi_integer,0,mpi_comm_world,code)

     if(rg==0)then
        ngll_mirror_tot = sum(sizes)
        displs(:) = sizes(:)
        do i = 2,Tdomain%n_proc
           displs(i) = displs(i)+displs(i-1)
        enddo
        do i = Tdomain%n_proc,2,-1
           displs(i) = displs(i-1)
        enddo
        displs(1) = 0
        allocate(glob_num_mirror_all(ngll_mirror_tot))
        allocate(x_gll_all(ngll_mirror_tot))
        allocate(y_gll_all(ngll_mirror_tot))
        allocate(z_gll_all(ngll_mirror_tot))
     endif

     call mpi_gatherv(x_gll,ngll_mirror,mpi_double_precision,x_gll_all,sizes,displs,mpi_double_precision,0,mpi_comm_world,code)
     call mpi_gatherv(y_gll,ngll_mirror,mpi_double_precision,y_gll_all,sizes,displs,mpi_double_precision,0,mpi_comm_world,code)
     call mpi_gatherv(z_gll,ngll_mirror,mpi_double_precision,z_gll_all,sizes,displs,mpi_double_precision,0,mpi_comm_world,code)
     call mpi_gatherv(glob_num_mirror,ngll_mirror,mpi_integer,glob_num_mirror_all,sizes,displs,mpi_integer,0,mpi_comm_world,code)

     if(rg==0)then
        allocate(index_num(ngll_mirror_tot))
        call indexx_i4b(glob_num_mirror_all,index_num)
        count_gll =1 
        do n = 2,ngll_mirror_tot
           i = index_num(n)
           if(glob_num_mirror_all(i)/=glob_num_mirror_all(index_num(n-1)))then
              count_gll = count_gll+1
           endif
        enddo
        open(48,file='mirror.vtk',status='unknown')    
        write(48,'(a)') '# vtk DataFile Version 2.0'
        write(48,'(a)') 'Mirror VTK file'
        write(48,'(a)') 'ASCII'
        write(48,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(48, '(a,i9,a)') 'POINTS ',count_gll, ' float'
        open(49,file='mirror.xyz',status='unknown')
        i = index_num(1)    
        write(48,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
        write(49,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
        do n = 2,ngll_mirror_tot
           i = index_num(n)
           if(glob_num_mirror_all(i)/=glob_num_mirror_all(index_num(n-1)))then
              write(48,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
              write(49,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
           endif
        enddo
        write(48,*)
        close(49)
        deallocate(sizes,displs,x_gll_all,y_gll_all,z_gll_all,glob_num_mirror_all,index_num)
     endif

     deallocate(glob_num_mirror)
     deallocate(x_gll)
     deallocate(y_gll)
     deallocate(z_gll)
     call mpi_barrier(mpi_comm_world,code)
     call mpi_finalize(code)
     write (*,*) "END ",rg
     stop
  endif ! end if output mirror glls

! ===================================================================================================

  allocate(glob_num_mirror(ngll_inner)) ! Check if the allocation is done properly
  allocate(x_gll(ngll_inner)) ! тот же самой
  allocate(y_gll(ngll_inner)) ! тот же самой
  allocate(z_gll(ngll_inner)) ! тот же самой
  ! 
  i = 0
  !
  call mpi_barrier(mpi_comm_world, code)

  !if (rg==0) then

  do n = 0,Tdomain%n_elem-1

     if(Tdomain%specel(n)%mirror_position == 1)then

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
                 ! if(Tdomain%specel(n)%win_mirror(x,y,z)/=0)then
                 i = i+1
                 iglob = Tdomain%specel(n)%Iglobnum(x,y,z)
                 glob_num_mirror(i) = (x+xs)+(y+ys)*nxtot+(z+zs)*(nxtot*nytot)
                 xa = Tdomain%Globcoord(0,iglob)
                 ya = Tdomain%Globcoord(1,iglob)
                 za = Tdomain%Globcoord(2,iglob)
                 x_gll(i) = Tdomain%rot(0,0)*xa + Tdomain%rot(0,1)*ya + Tdomain%rot(0,2)*za       
                 y_gll(i) = Tdomain%rot(1,0)*xa + Tdomain%rot(1,1)*ya + Tdomain%rot(1,2)*za
                 z_gll(i) = Tdomain%rot(2,0)*xa + Tdomain%rot(2,1)*ya + Tdomain%rot(2,2)*za
                 ! endif
              enddo
           enddo
        enddo

     endif

  enddo

  call mpi_barrier(mpi_comm_world, code)

  if(rg==0)then
     allocate(sizes(Tdomain%n_proc))
     allocate(displs(Tdomain%n_proc))
  endif

  call mpi_barrier(mpi_comm_world, code)

  call mpi_gather(ngll_inner,1,mpi_integer,sizes,1,mpi_integer,0,mpi_comm_world,code)

  if(rg==0)then
     ngll_mirror_tot = sum(sizes)
     displs(:) = sizes(:)
     do i = 2,Tdomain%n_proc
        displs(i) = displs(i)+displs(i-1)
     enddo
     do i = Tdomain%n_proc,2,-1
        displs(i) = displs(i-1)
     enddo
     displs(1) = 0
     allocate(glob_num_mirror_all(ngll_mirror_tot))
     allocate(x_gll_all(ngll_mirror_tot))
     allocate(y_gll_all(ngll_mirror_tot))
     allocate(z_gll_all(ngll_mirror_tot))
  endif

  call mpi_gatherv(x_gll,ngll_inner,mpi_double_precision,x_gll_all,sizes,displs,mpi_double_precision,0,mpi_comm_world,code)
  call mpi_gatherv(y_gll,ngll_inner,mpi_double_precision,y_gll_all,sizes,displs,mpi_double_precision,0,mpi_comm_world,code)
  call mpi_gatherv(z_gll,ngll_inner,mpi_double_precision,z_gll_all,sizes,displs,mpi_double_precision,0,mpi_comm_world,code)
  call mpi_gatherv(glob_num_mirror,ngll_inner,mpi_integer,glob_num_mirror_all,sizes,displs,mpi_integer,0,mpi_comm_world,code)
  call mpi_barrier(mpi_comm_world, code)

  ! <SA> Second sorting and scattering the line in the specfem mirror file to every proc
  ! so that merging won't be needed afterwards for the convolution

  allocate(Tdomain%index_pos(ngll_inner))

  if(rg==0) then
     allocate(index_num(ngll_mirror_tot))
     allocate(index_pos_all(ngll_mirror_tot))
     call indexx_i4b(glob_num_mirror_all,index_num)
     ! call indexx_i4b(index_num, index_pos_all)
     count_gll = 1
     index_pos_all(1) = count_gll
     do n = 2,ngll_mirror_tot
        i = index_num(n)
        if(glob_num_mirror_all(i)/=glob_num_mirror_all(index_num(n-1)))then
           count_gll = count_gll+1
        endif
        index_pos_all(i) = count_gll
     enddo
  endif
  
  call mpi_scatterv(index_pos_all,sizes,displs,mpi_integer,Tdomain%index_pos,ngll_inner,mpi_integer,0,mpi_comm_world,code)

  if(Tdomain%t_reversal_mirror==5)then ! Output inner mirror's GLLs for globalwave propagation solver 
     if(rg==0)then
        print*, 'Entered the third if for root'
        ! <SA> next lines computed before so the double sorting can be computed
        ! allocate(index_num(ngll_mirror_tot))
        ! call indexx_i4b(glob_num_mirror_all,index_num)
        count_gll =1 
        do n = 2,ngll_mirror_tot
           i = index_num(n)
           if(glob_num_mirror_all(i)/=glob_num_mirror_all(index_num(n-1)))then
              count_gll = count_gll+1
           endif
        enddo
        open(48,file='mirror.vtk',status='unknown')    
        write(48,'(a)') '# vtk DataFile Version 2.0'
        write(48,'(a)') 'Mirror VTK file'
        write(48,'(a)') 'ASCII'
        write(48,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(48, '(a,i9,a)') 'POINTS ',count_gll, ' float'
        open(49,file='mirror.xyz',status='unknown')
        i = index_num(1)
        ii = minval(glob_num_mirror_all)
        write(48,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
        write(49,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
        do n = 2,ngll_mirror_tot
           i = index_num(n)
           if(glob_num_mirror_all(i)/=glob_num_mirror_all(index_num(n-1)))then
              write(48,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
              write(49,*)real(x_gll_all(i)),real(y_gll_all(i)),real(z_gll_all(i))
           endif
        enddo
        write(48,*)
        close(49)
        deallocate(sizes,displs,x_gll_all,y_gll_all,z_gll_all,glob_num_mirror_all,index_num)
     endif
     
     deallocate(glob_num_mirror)
     deallocate(x_gll)
     deallocate(y_gll)
     deallocate(z_gll)
     call mpi_barrier(mpi_comm_world,code)
     call mpi_finalize(code)
     write (*,*) "END ",rg
     stop
  endif ! end if output inner mirror glls

  return
end subroutine mirror_definition

