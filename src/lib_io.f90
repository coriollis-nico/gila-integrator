module lib_io
  !! Utilities for file output
  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128
  implicit none

private
public matrix_to_file
public safe_open
public save_integral_conditions

interface matrix_to_file
  !! For writing 2d matrices to files
  module procedure :: realrow_to_file
  module procedure :: realmatrix_to_file
  module procedure :: intmatrix_to_file
end interface matrix_to_file

contains

subroutine realrow_to_file(matrix, file_id)

  real(qp), intent(in), dimension(:) :: matrix
  integer, intent(in) :: file_id

  integer :: j

  do j = 1, size(matrix, 1)
    write(file_id, *) matrix(j)
  end do

end subroutine realrow_to_file

subroutine realmatrix_to_file(matrix, file_id)

  real(qp), intent(in), dimension(:,:) :: matrix
  integer, intent(in) :: file_id

  integer :: j

  do j = 1, size(matrix, 1)
    write(file_id, *) matrix(j, :)
  end do

end subroutine realmatrix_to_file

subroutine intmatrix_to_file(matrix, file_id)

  integer, intent(in), dimension(:,:) :: matrix
  integer, intent(in) :: file_id

  integer :: j

  do j = 1, size(matrix, 1)
    write(file_id, *) matrix(j, :)
  end do

end subroutine intmatrix_to_file

subroutine safe_open(file_name, file_id, file_dir)

  !! Ensures file is not overwriten when opened

  character(len=*), intent(in) :: file_name
  !! Name of output file
  integer, intent(inout) :: file_id
  !! File id.
  character(len=*), intent(in), optional :: file_dir
  !! Directory of output file

  logical :: exists


    select case( present(file_dir) )
      case(.true.)
        inquire(file=file_dir//"/"//file_name, exist=exists)
        if (.not. exists) then
          call execute_command_line("mkdir --parent "//file_dir, wait=.true.)
          open(newunit=file_id, file=file_dir//"/"//file_name, &
                & status="new", action="write")
        else
          stop "Preventing overwrite"
        end if
      case(.false.)
        inquire(file=file_name, exist=exists)
        if (.not. exists) then
          open(newunit=file_id, file=file_name, status="new", action="write")
        else
          stop "Preventing overwrite"
        end if
    end select


end subroutine safe_open

subroutine save_integral_conditions(x, y0, file_id)

  !! Saves numeric integration parameters to file

  real(qp), intent(in) :: x(:)
  real(qp), intent(in) :: y0
  integer, intent(inout) :: file_id

  write(file_id, *) "# ---"
  write(file_id, *) "# x:", x(1), "->", x(size(x))
  write(file_id, *) "# y0:", y0
  write(file_id, *) "# n = ", size(x)
  write(file_id, *) "# h = ", x(2) - x(1)

end subroutine save_integral_conditions

end module lib_io
