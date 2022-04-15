function blocksize(bs) bind(C, name='blocksize')
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: bs
  integer(C_INT) :: val
  integer(C_INT) :: blocksize
  integer, external :: ilaenv
  
  val = ilaenv(1, 'DGEQRF', '', bs, bs, bs, bs)
  
  blocksize = val
end function blocksize
