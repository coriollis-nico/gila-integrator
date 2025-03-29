program curve_sandwich
  !! For specified [[gila_conditions:m]], [[gila_conditions:p]], [[gila_conditions:p]]
  !! combinations calculates \( y \) and \( \epsilon_{1 \to 4} \) slow-roll parameters.
  !! Also limit curves.
  !!
  !! Also calculates slow-roll curves for each solution
  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128,&
          stdout => output_unit
  use gila
  use lib_io
  implicit none

  integer :: i, k, j

  integer, dimension(4), parameter  :: mt = [ 3,         8,         10,       99 ]
  !! [[gila_conditions:m]]
  integer, dimension(4), parameter  :: pt = [ 1,         2,         10,       99 ]
  !! [[gila_conditions:p]]
  real(qp), dimension(3), parameter :: lt = [ 1.e-17_qp, 1.e-22_qp, 1.e-27_qp ]
  !! [[gila_conditions:l]]

  real(qp), parameter :: xi = 0.0_qp
  !! Initial [[gila_friedmann:x]] value
  real(qp), parameter :: xf = -100.0_qp
  !! Final [[gila_friedmann:y]] value
  real(qp), parameter :: yi = 1.0_qp
  !! Initial [[gila_friedmann:y]] value
  integer, parameter :: n = 20000
  !! Number of integration steps

  real(qp), dimension(n), parameter :: x = [( i * (xf - xi)/n, i = 0, n-1 )]
  !! \( x \) values array
  real(qp), dimension(n) :: y
  !! Solution array
  real(qp), dimension(n, 5) :: slowroll_data
  !! \( \epsilon_0 \) array

  type(gila_conditions) :: conditions
  real(qp), parameter :: matter_density = 0.31_qp
  real(qp), parameter :: rad_density = 8.4e-5_qp
  real(qp), parameter :: dark_density = 0.69_qp

  character(len=*), parameter :: data_dir = "data/sims/curve_sandwich"
  character(len=32) :: filename

  integer :: x_grid
  integer :: out_id
  integer :: mp_id, l_id, lim_id
  integer :: sr
! ------------------------------------------------------------------------------------ !

call safe_open("mp.dat", mp_id, file_dir=data_dir)

  call save_integral_conditions(x, yi, mp_id)
  do i = 1, 4
    write(mp_id, *) mt(i), pt(i)
  end do

close(mp_id)


call safe_open("l.dat", l_id, file_dir=data_dir)

  call save_integral_conditions(x, yi, mp_id)
  do i = 1, 3
    write(mp_id, *) lt(i)
  end do

close(l_id)


! --- x


call safe_open("x_grid.dat", x_grid, file_dir=data_dir)

  call save_integral_conditions(x, yi, x_grid)

  do k = 1, size(x)
    write(x_grid, *) x(k)
  end do

close(x_grid)

! ---

conditions%cosmos="early"
conditions%beta=0.0_qp
conditions%l_tilde=0.0_qp
conditions%r=0
conditions%s=0
conditions%Omega_M = matter_density
conditions%Omega_R = rad_density
conditions%Omega_dark = dark_density


do k = 1, 3

  conditions%l=lt(k)

  do i = 1, 4

  conditions%m=mt(i)
  conditions%p=pt(i)

  write(filename, '(2(a1, i2.2), a1, es7.1e2)') 'm', conditions%m, &
                                                'p', conditions%p, &
                                                'l', conditions%l

  call safe_open(trim(filename)//".dat", out_id, file_dir=data_dir)

    y = gila_solution(x, yi, conditions)

    call save_gila_genconditions(conditions, out_id)
    call save_gila_genconditions(conditions, stdout)
    call save_integral_conditions(x, yi, out_id)

    call matrix_to_file(y, out_id)

  close(out_id)

  !v - slowroll
  write(filename, '(a5, 2(a1, i2.2), a1, es7.1e2)') 'slow-', 'm', conditions%m, &
                                                'p', conditions%p, &
                                                'l', conditions%l

  call safe_open(trim(filename)//".dat", sr, file_dir=data_dir)

    slowroll_data(:,1:2) = slowroll0(x, y, conditions%l)
    do j = 3, 5
      slowroll_data(:,j) = slowroll(slowroll_data(:,1), slowroll_data(:,j-1))
    end do

    write(sr, *) "# N   ϵ0   ϵ1   ϵ2   ϵ3"
    call matrix_to_file(slowroll_data, sr)

  close(sr)

  end do

end do

! limit curves ----------------------------------------------------------------------------------- !

do k = 1, 3

  conditions%cosmos="early"
  conditions%l=lt(k)

  write(filename, '(a5, es7.1e2)') 'lim_l', conditions%l

  call safe_open(trim(filename)//".dat", lim_id, file_dir=data_dir)

    y = gila_solution_limit(x, yi, conditions)

    call save_gila_limconditions(conditions, lim_id)
    call save_gila_limconditions(conditions, stdout)
    call save_integral_conditions(x, yi, lim_id)

    call matrix_to_file(y, lim_id)

  close(lim_id)

end do

end program curve_sandwich
