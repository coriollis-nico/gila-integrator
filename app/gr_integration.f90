program gr_integration
  !! Simulates a GR universe from initial conditions
  use iso_fortran_env,  &
    only: qp => real128
  use gr
  use io
  implicit none

! integration conditions
  real(qp), parameter :: x_i = 0.0_qp
!! Initial [[gila_friedmann:x]] value
  real(qp), parameter :: x_f = -50.0_qp
!! Final [[gila_friedmann:x]] value
  real(qp), parameter :: yi = 1.0_qp
!! Initial [[gila_friedmann:y]] value
  integer, parameter :: n = 10000
!! Number of integration steps
  type(gr_conditions), parameter :: cond = gr_conditions_default
!! Specify GR conditions

  character(*), parameter :: data_dir = "data/sims/gr"
  character(*), parameter :: data_file = "gr.dat"
!! Where to save data (relative to project root dir)
  integer                 :: io_gr
!! `data_dir` identifier

  integer :: i
!! Iterator

  real(qp), dimension(n), parameter :: x = [( x_i + i*(x_f - x_i)/(n), i = 0, n-1 )]
!! \( x \) values array
  real(qp), dimension(n) :: y
!! solution

! ------------------------------------------------------------------------------------ !

  print '(a)', "Calculating solution..."
  y = gr_solution(x, yi, cond)

  print '(a)', "Saving data..."
  call safe_open(data_file, io_gr, data_dir)
  write(io_gr, '(a)') "# x  y"
  do i = 1, n
    write(io_gr, *) x(i), y(i)
  end do
  close(io_gr)

end program gr_integration
