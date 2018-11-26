module s_src_ext

type :: src_ext

integer :: n_gll
integer :: spln_order
integer :: nt
integer :: ngll_mirror,ngll_mirror_tot
integer :: ind_old,offset
real :: dt,t0
logical :: init
!real, dimension(:), pointer :: 
integer, dimension(:), pointer :: displs,sizes,glob_num_mirror,index_num,glob_num_mirror_all
real, dimension(:), pointer :: tmp_x,tmp_x_all,tmp_y,tmp_y_all,tmp_z,tmp_z_all
real, dimension(:,:,:), pointer :: b_coefs
end type

end module
