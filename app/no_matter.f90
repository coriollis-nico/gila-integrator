program no_matter
  !! For specified [[gila_conditions:m]], [[gila_conditions:p]]
  !! combinations calculates same evolution as [[curve_sandwich]] for conditions starting on
  !! inflation exit

  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128,&
          stdout => output_unit
  use gila
  use lib_io
  implicit none

  integer :: i, k, j

  integer, dimension(3), parameter  :: mt = [ 3,         8,         10 ]
  !! [[gila_conditions:m]]
  integer, dimension(3), parameter  :: pt = [ 1,         2,         10 ]
  !! [[gila_conditions:p]]

  real(qp), parameter :: xi = 0.0_qp
  !! Initial \( \bar{a}_{(i)}} \) value
  real(qp), parameter :: xf = -60.0_qp
  !! Final \( \bar{a}_{(i)}} \) value
  real(qp), parameter :: yi = 1.0_qp
  !! Initial \( \bar{H}_{(i)} \) value
  integer, parameter :: n = 20000
  !! Number of integration steps

  real(qp), dimension(n), parameter :: x = [( i * (xf - xi)/n, i = 0, n-1 )]
  !! ´x´ values array
  real(qp), dimension(n) :: y
  !! Solution array
  real(qp), dimension(n, 5) :: slowroll_data
  !! \( \epsilon_n \) array

  type(gila_conditions) :: conditions
  real(qp), parameter :: matter_density = 0.0_qp
  real(qp), parameter :: rad_density = 0.31_qp
  real(qp), parameter :: dark_density = 0.69_qp

  character(len=*), parameter :: data_dir = "data/sims/no_matter"
  character(len=32) :: filename

  integer :: x_grid
  integer :: out_id
  integer :: mp_id, lim_id
  integer :: sr
! ------------------------------------------------------------------------------------ !

call safe_open("mp.dat", mp_id, file_dir=data_dir)

  call save_integral_conditions(x, yi, mp_id)
  do i = 1, 3
    write(mp_id, *) mt(i), pt(i)
  end do

close(mp_id)

! --- x


call safe_open("x_grid.dat", x_grid, file_dir=data_dir)

  call save_integral_conditions(x, yi, x_grid)

  do k = 1, size(x)
    write(x_grid, *) x(k)
  end do

close(x_grid)

! ---

conditions%cosmos="einf"
conditions%beta=0.0_qp
conditions%l=1.0_qp
conditions%l_tilde=0.0_qp
conditions%r=0
conditions%s=0
conditions%Omega_M = matter_density
conditions%Omega_R = rad_density
conditions%Omega_dark = dark_density


do i = 1, 3

  conditions%m=mt(i)
  conditions%p=pt(i)

  write(filename, '(2(a1, i2.2))') 'm', conditions%m, &
                                  'p', conditions%p

  call safe_open(trim(filename)//".dat", out_id, file_dir=data_dir)

  y = gila_solution(x, yi, conditions)

  call save_gila_genconditions(conditions, out_id)
  call save_gila_genconditions(conditions, stdout)
  call save_integral_conditions(x, yi, out_id)

  call matrix_to_file(y, out_id)

  close(out_id)

  !v - slowroll
  write(filename, '(a5, 2(a1, i2.2))') 'slow-', 'm', conditions%m, &
                                                'p', conditions%p

  call safe_open(trim(filename)//".dat", sr, file_dir=data_dir)

  slowroll_data(:,1) = x
  slowroll_data(:,2) = 1.0_qp/y
  do j = 3, 5
    slowroll_data(:,j) = slowroll(slowroll_data(:,1), slowroll_data(:,j-1))
  end do

  write(sr, *) "# N   ϵ0   ϵ1   ϵ2   ϵ3"
  call matrix_to_file(slowroll_data, sr)

  close(sr)

end do

end program no_matter
