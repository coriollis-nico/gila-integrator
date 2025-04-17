program gr_integration
  !! Simulates a GR universe from initial conditions
  use iso_fortran_env,  &
    only: sp => real32, &
    dp => real64, &
    qp => real128
  use gila
  use lib_io
  implicit none

! integration conditions
  real(qp), parameter :: x_i = 0.0_qp
!! Initial [[gila_friedmann:x]] value
  real(qp), parameter :: x_f = -50.0_qp
!! Final [[gila_friedmann:y]] value
  real(qp), parameter :: yi = 1.0_qp
!! Initial [[gila_friedmann:y]] value
  integer, parameter :: n = 10000
!! Number of integration steps
  type(gila_conditions), parameter :: gr_conditions = gila_grdefault
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

  print *, "Calculating solution..."
  y = gila_solution(x, yi, gr_conditions)

  print *, "Saving data"
  call safe_open(data_file, io_gr, data_dir)
  call save_gila_genconditions(gr_conditions, io_gr)
  call save_integral_conditions(x, yi, io_gr)
  write(io_gr, *) "# x  y"
  do i = 1, n
    write(io_gr, *) x(i), y(i)
  end do
  close(io_gr)

end program gr_integration
