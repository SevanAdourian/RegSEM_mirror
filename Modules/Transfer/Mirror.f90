module smirror

type :: mirror

integer*8 :: recl_mirror
integer :: decim_fact
integer :: spln_order
integer :: count_rec
integer :: lunit
integer :: lunit2
integer :: ngll_mirror_tot
doubleprecision, dimension(:,:), pointer :: tmp
doubleprecision, dimension(:), pointer :: b

end type 

end module smirror
