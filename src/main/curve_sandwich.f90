program curve_sandwich
  !! For specified [[gila_conditions:m]], [[gila_conditions:p]] combinations calculates
  !! \( y \) and \( \epsilon_{1 \to 4} \) slow-roll parameters
  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128,&
          stdout => output_unit
  use gila
  use lib_io
  implicit none

  integer :: i, k

  integer, dimension(3), parameter  :: m = [ 3,         8,         10 ]
  !! [[gila_conditions:m]]
  integer, dimension(3), parameter  :: p = [ 1,         2,         10 ]
  !! [[gila_conditions:p]]
  real(qp), dimension(3), parameter :: l = [ 1.e-17_qp, 1.e-22_qp, 1.e-27_qp ]
  !! [[gila_conditions:l]]

  real(qp), parameter :: xi = 0.0_qp
  !! Initial [[gila_friedmann:x]] value
  real(qp), parameter :: xf = -100.0_qp
  !! Final [[gila_friedmann:y]] value
  real(qp), parameter :: yi = 0.0_qp
  !! Initial [[gila_friedmann:y]] value
  integer, parameter :: n = 5000
  !! Number of integration steps
  type(gila_conditions) :: conditions
  !! Specify GILA conditions. Specific [[gila_conditions:m]], [[gila_conditions:p]],
  !! [[gila_conditions:l]] specified per iteration

  real(qp), dimension(n), parameter :: x = [( i * (xf - xi)/n, i = 0, n-1 )]
  !! \( x \) values array
  real(qp), dimension(n, 5) :: y
  !! - `y(1,:)` stores the solution
  !! - `y(j,:)` stores the j-1th slowroll param.

  character(len=*), parameter :: data_dir = "data/sims/curve_sandwich"

  integer :: x_grid
  integer :: m03p01l17, m08p02l17, m10p10l17
  integer :: m03p01l22, m08p02l22, m10p10l22
  integer :: m03p01l27, m08p02l27, m10p10l27

! ------------------------------------------------------------------------------------ !

conditions%cosmos = "early"
conditions%beta = 0.0_qp
conditions%l_tilde = 0.0_qp

! --- x

call safe_open("x_grid.dat", x_grid, file_dir=data_dir)

  call save_integral_conditions(x, yi, x_grid)

  do k = 1, size(x)
    write(x_grid, *) x(k)
  end do

close(x_grid)

! --- m(1), p(1), l(...)

call safe_open("m03p01l17.dat", m03p01l17, file_dir=data_dir)

  conditions%m = m(1)
  conditions%p = p(1)
  conditions%l = l(1)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m03p01l17)
  call save_integral_conditions(x, yi, m03p01l17)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m03p01l17)

close(m03p01l17)

call safe_open("m03p01l22.dat", m03p01l22, file_dir=data_dir)

  conditions%m = m(1)
  conditions%p = p(1)
  conditions%l = l(2)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m03p01l22)
  call save_integral_conditions(x, yi, m03p01l22)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m03p01l22)

close(m03p01l22)

call safe_open("m03p01l27.dat", m03p01l27, file_dir=data_dir)

  conditions%m = m(1)
  conditions%p = p(1)
  conditions%l = l(3)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m03p01l27)
  call save_integral_conditions(x, yi, m03p01l27)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m03p01l27)

close(m03p01l27)


! --- m(2), p(2), l(...)

call safe_open("m08p02l17.dat", m08p02l17, file_dir=data_dir)

  conditions%m = m(2)
  conditions%p = p(2)
  conditions%l = l(1)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m08p02l17)
  call save_integral_conditions(x, yi, m08p02l17)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m08p02l17)

close(m08p02l17)


call safe_open("m08p02l22.dat", m08p02l22, file_dir=data_dir)

  conditions%m = m(2)
  conditions%p = p(2)
  conditions%l = l(2)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m08p02l22)
  call save_integral_conditions(x, yi, m08p02l22)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m08p02l22)

close(m08p02l22)

call safe_open("m08p02l27.dat", m08p02l27, file_dir=data_dir)

  conditions%m = m(2)
  conditions%p = p(2)
  conditions%l = l(3)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m08p02l27)
  call save_integral_conditions(x, yi, m08p02l27)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m08p02l27)

close(m08p02l27)


! --- m(3), p(4), l(...)

call safe_open("m10p10l17.dat", m10p10l17, file_dir=data_dir)

  conditions%m = m(3)
  conditions%p = p(3)
  conditions%l = l(1)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m10p10l17)
  call save_integral_conditions(x, yi, m10p10l17)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m10p10l17)

close(m10p10l17)


call safe_open("m10p10l22.dat", m10p10l22, file_dir=data_dir)

  conditions%m = m(3)
  conditions%p = p(3)
  conditions%l = l(2)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m10p10l22)
  call save_integral_conditions(x, yi, m10p10l22)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m10p10l22)

close(m10p10l22)

call safe_open("m10p10l27.dat", m10p10l27, file_dir=data_dir)

  conditions%m = m(3)
  conditions%p = p(3)
  conditions%l = l(3)

  call save_gila_genconditions(conditions, stdout)

  call save_gila_genconditions(conditions, m10p10l27)
  call save_integral_conditions(x, yi, m10p10l27)

  y = gila_solution(x, yi, 4, conditions)

  call matrix_to_file(y, m10p10l27)

close(m10p10l27)


end program curve_sandwich
