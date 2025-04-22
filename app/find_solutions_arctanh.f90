program find_solutions_arctanh
  !! TODO Write
  use iso_fortran_env,  &
    only: sp => real32, &
    dp => real64, &
    qp => real128
  use gi_arctanh
  use lib_io
  implicit none

  integer :: i, k

  real(qp), dimension(3), parameter :: lt = [ 1.e-17_qp, 1.e-22_qp, 1.e-27_qp ]
  !! [[gi_actanh_conditions:l]]

  real(qp), parameter :: xi = 0.0_qp
  !! Initial [[gi_arctanh_friedmann:x]] value
  real(qp), parameter :: xf = -100.0_qp
  !! Final [[gi_arctanh_friedmann:y]] value
  real(qp), parameter :: yi = 1.0_qp
  !! Initial [[gi_arctanh_friedmann:y]] value
  integer, parameter :: n = 20000
  !! Number of integration steps

  real(qp), dimension(n), parameter :: x = [( i * (xf - xi)/n, i = 0, n-1 )]
  !! \( x \) values array
  real(qp), dimension(n, size(lt)) :: y
  !! Solutions array
  real(qp), dimension(n, size(lt), 2) :: sr_n_0
  !! \( N, \epsilon_0 \) array
  real(qp), dimension(n, size(lt)) :: sr1, sr2, sr3
  !! \( \epsilon \) array

  type(gi_arctanh_conditions) :: conditions
  real(qp), parameter :: matter_density = 0.31_qp
  real(qp), parameter :: rad_density = 8.4e-5_qp
  real(qp), parameter :: dark_density = 0.69_qp

  character(len=*), parameter :: data_dir = "data/sims/solutions_arctanh"

  integer :: x_id
  integer :: y_id
  integer :: l_id
  integer :: srn_id, sr0_id, sr1_id, sr2_id, sr3_id
! ------------------------------------------------------------------------------------ !

  conditions%Omega_M = matter_density
  conditions%Omega_R = rad_density
  conditions%Omega_dark = dark_density

  call safe_open("l.dat", l_id, file_dir=data_dir)

  do k = 1, size(lt)
    write(l_id, *) lt(k)
  end do

  close(l_id)


! --- x

  call safe_open("x.dat", x_id, file_dir=data_dir)

  call matrix_to_file(x, x_id)

  close(x_id)


! ---

  do k = 1, 3

    conditions%l=lt(k)

    y(:, k) = gi_arctanh_solution(x, yi, conditions)

    sr_n_0(:, k, :) = gi_arctanh_sr0(x, y(:, k), conditions%l)
    sr1(:, k) = gi_arctanh_sr(sr_n_0(:, k, 1), sr_n_0(:, k, 2))
    sr2(:, k) = gi_arctanh_sr(sr_n_0(:, k, 1), sr1(:, k))
    sr3(:, k) = gi_arctanh_sr(sr_n_0(:, k, 1), sr2(:, k))

  end do

  call safe_open("y.dat", y_id, file_dir=data_dir)
  do k = 1, size(lt)
    write(y_id, *) y(:, k)
  end do
  close(y_id)

  call safe_open("n.dat", srn_id, file_dir=data_dir)
  do k = 1, size(lt)
    write(srn_id, *) sr_n_0(:, k, 1)
  end do
  close(srn_id)

  call safe_open("sr0.dat", sr0_id, file_dir=data_dir)
  do k = 1, size(lt)
    write(sr0_id, *) sr_n_0(:, k, 2)
  end do
  close(sr0_id)

  call safe_open("sr1.dat", sr1_id, file_dir=data_dir)
  do k = 1, size(lt)
    write(sr1_id, *) sr1(:, k)
  end do
  close(sr1_id)

  call safe_open("sr2.dat", sr2_id, file_dir=data_dir)
  do k = 1, size(lt)
    write(sr2_id, *) sr2(:, k)
  end do
  close(sr2_id)

  call safe_open("sr3.dat", sr3_id, file_dir=data_dir)
  do k = 1, size(lt)
    write(sr3_id, *) sr3(:, k)
  end do
  close(sr3_id)


end program find_solutions_arctanh
