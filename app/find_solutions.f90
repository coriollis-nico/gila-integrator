program find_solutions
  !! For specified [[gila_conditions:m]], [[gila_conditions:p]], [[gila_conditions:l]]
  !! combinations calculates \( y \) and \( \epsilon_{1 \to 4} \) slow-roll parameters.
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

  integer :: i, k, j, x_c, e_c

  integer, dimension(3), parameter  :: mt = [ 94, 95, 96 ]
  !! [[gila_conditions:m]]
  integer, dimension(3), parameter  :: pt = [ 3, 4, 5 ]
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
  real(qp), dimension(n, size(mt), size(pt), size(lt)) :: y
  !! Solutions array
  real(qp), dimension(n, size(mt), size(pt), size(lt), 2) :: sr_n_0
  !! \( N, \epsilon_0 \) array
  real(qp), dimension(n, size(mt), size(pt), size(lt)) :: sr1, sr2, sr3
  !! \( \epsilon \) array

  type(gila_conditions) :: conditions
  real(qp), parameter :: matter_density = 0.31_qp
  real(qp), parameter :: rad_density = 8.4e-5_qp
  real(qp), parameter :: dark_density = 0.69_qp

  character(len=*), parameter :: data_dir = "data/sims/extra_cases"
  character(len=32) :: filename

  integer :: x_id
  integer :: y_id
  integer :: mpl_id
  integer :: srn_id, sr0_id, sr1_id, sr2_id, sr3_id
! ------------------------------------------------------------------------------------ !

call safe_open("mpl.dat", mpl_id, file_dir=data_dir)

  do i = 1, size(mt)
  do j = 1, size(pt)
  do k = 1, size(lt)
    write(mpl_id, *) mt(i), pt(j), lt(k)
  end do
  end do
  end do

close(mpl_id)


! --- x

call safe_open("x.dat", x_id, file_dir=data_dir)

  call matrix_to_file(x, x_id)

close(x_id)


! ---

conditions%cosmos="early"
conditions%beta=0.0_qp
conditions%l_tilde=0.0_qp
conditions%r=0
conditions%s=0
conditions%Omega_M = matter_density
conditions%Omega_R = rad_density
conditions%Omega_dark = dark_density


do i = 1, size(mt)

  conditions%m=mt(i)

  do j = 1, size(pt)

    conditions%p=pt(j)

    do k = 1, size(lt)

      conditions%l=lt(k)

      y(:, i, j, k) = gila_solution(x, yi, conditions)

      sr_n_0(:, i, j, k, :) = slowroll0(x, y(:, i, j, k), conditions%l)
      sr1(:, i, j, k) = slowroll(sr_n_0(:, i, j, k, 1), sr_n_0(:, i, j, k, 2))
      sr2(:, i, j, k) = slowroll(sr_n_0(:, i, j, k, 1), sr1(:, i, j, k))
      sr3(:, i, j, k) = slowroll(sr_n_0(:, i, j, k, 1), sr2(:, i, j, k))

    end do

  end do

end do

call safe_open("y.dat", y_id, file_dir=data_dir)
  do i = 1, size(mt)
  do j = 1, size(pt)
  do k = 1, size(lt)
    write(y_id, *) y(:, i, j, k)
  end do
  end do
  end do
close(y_id)

call safe_open("n.dat", srn_id, file_dir=data_dir)
  do i = 1, size(mt)
  do j = 1, size(pt)
  do k = 1, size(lt)
    write(srn_id, *) sr_n_0(:, i, j, k, 1)
  end do
  end do
  end do
close(srn_id)

call safe_open("sr1.dat", sr1_id, file_dir=data_dir)
  do i = 1, size(mt)
  do j = 1, size(pt)
  do k = 1, size(lt)
    write(sr1_id, *) sr1(:, i, j, k)
  end do
  end do
  end do
close(sr1_id)

call safe_open("sr2.dat", sr2_id, file_dir=data_dir)
  do i = 1, size(mt)
  do j = 1, size(pt)
  do k = 1, size(lt)
    write(sr2_id, *) sr2(:, i, j, k)
  end do
  end do
  end do
close(sr2_id)

call safe_open("sr3.dat", sr3_id, file_dir=data_dir)
  do i = 1, size(mt)
  do j = 1, size(pt)
  do k = 1, size(lt)
    write(sr3_id, *) sr3(:, i, j, k)
  end do
  end do
  end do
close(sr3_id)


end program find_solutions
