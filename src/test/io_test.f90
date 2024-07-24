program io_test
  !! Test for [[lib_io]] module
  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128,&
          io => output_unit
  use lib_io
  implicit none

  integer,  dimension(4,4) :: mi
  real(qp), dimension(4,4) :: mr

  integer :: j,k

  character(*), parameter :: separator = &
    &"-----------------------------------------&
    &-----------------------------------------"

! ------------------------------------------------------------------------------------ !

do j = 1, 4
  do k = 1, 4
    mi(j,k) = j + k
    mr(j,k) = exp(1.0_qp * j * k)
  end do
end do

! matrix_to_file

call matrix_to_file(mi, io)
print *, separator
call matrix_to_file(mr, io)

end program io_test
