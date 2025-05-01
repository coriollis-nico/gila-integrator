program find_solutions_initial
  !! For specified [[gi_conditions:m]], [[gi_conditions:p]], [[gi_conditions:l]]
  !! combinations calculates \( y \) and \( \epsilon_{1 \to 4} \) slow-roll parameters.
  !!
  !! Also calculates slow-roll curves for each solution
  use iso_fortran_env,   &
    only: qp => real128, &
    stdout => output_unit
  use gi
  use io
  implicit none

  integer :: i, k

  integer, dimension(3), parameter  :: mt = [ 3, 8, 10 ]
  !! [[gi_conditions:m]]
  integer, dimension(size(mt)), parameter  :: pt = [ 1, 2, 10 ]
  !! [[gi_conditions:p]]
  real(qp), dimension(3), parameter :: lt = [ 1.e-17_qp, 1.e-22_qp, 1.e-27_qp ]
  !! [[gi_conditions:l]]

  real(qp), parameter :: xi = 0.0_qp
  !! Initial [[gi_friedmann:x]] value
  real(qp), parameter :: xf = -100.0_qp
  !! Final [[gi_friedmann:y]] value
  real(qp), parameter :: yi = 1.0_qp
  !! Initial [[gi_friedmann:y]] value
  integer, parameter :: n = 20000
  !! Number of integration steps

  real(qp), dimension(n), parameter :: x = [( i * (xf - xi)/n, i = 0, n-1 )]
  !! \( x \) values array
  real(qp), dimension(n, size(mt), size(lt)) :: y
  !! Solutions array
  real(qp), dimension(n, size(mt), size(lt), 2) :: sr_n_0
  !! \( N, \epsilon_0 \) array
  real(qp), dimension(n, size(mt), size(lt)) :: sr1, sr2, sr3
  !! \( \epsilon \) array

  real(qp), parameter :: matter_density = 0.31_qp
  real(qp), parameter :: rad_density = 8.4e-5_qp
  real(qp), parameter :: dark_density = 0.69_qp

  character(len=*), parameter :: data_dir = "data/sims/solutions_initial"

  integer :: x_id
  integer :: y_id
  integer :: mpl_id
  integer :: srn_id, sr1_id, sr2_id, sr3_id
! ------------------------------------------------------------------------------------ !

  write(io_str, '((a), i5, (a))') '(', n, '(es30.20))'

  print '(a)', io_str

  print '(a)', "Saving mpl data..."
  call safe_open("mpl.dat", mpl_id, file_dir=data_dir)

  do k = 1, size(lt)
    do i = 1, size(mt)
      write(mpl_id, '(3(i4))') mt(i), pt(i), int(log10(lt(k)))
    end do
  end do

  close(mpl_id)


! --- x
    print '(a)', "Saving x data..."
  call safe_open("x.dat", x_id, file_dir=data_dir)
  call matrix_to_file(x, x_id)
  close(x_id)


! ---

  print '(a)', "Finding solutions"
  do concurrent (i = 1:size(mt))
    do concurrent (k = 1:size(lt))
      y(:, i, k) = gi_solution(x, yi, gi_conditions(&
                                        ld = 1.0_qp, &
                                        l = lt(k), &
                                        p = pt(i), &
                                        m = mt(i), &
                                        Omega_M = matter_density, &
                                        Omega_R = rad_density, &
                                        Omega_dark = dark_density))
    end do
  end do

  print '(a)', "Finding ϵ data..."
  do concurrent (i = 1:size(mt))
    do concurrent (k = 1:size(lt))
      sr_n_0(:, i, k, :) = gi_sr0(x, y(:, i, k), gi_conditions(&
                                                  ld = 1.0_qp, &
                                                  l = lt(k), &
                                                  p = pt(i), &
                                                  m = mt(i), &
                                                  Omega_M = matter_density, &
                                                  Omega_R = rad_density, &
                                                  Omega_dark = dark_density))
      sr1(:, i, k) = gi_sr(sr_n_0(:, i, k, 1), sr_n_0(:, i, k, 2))
      sr2(:, i, k) = gi_sr(sr_n_0(:, i, k, 1), sr1(:, i, k))
      sr3(:, i, k) = gi_sr(sr_n_0(:, i, k, 1), sr2(:, i, k))
    end do
  end do


    print '(a)', "Saving data..."
  call safe_open("y.dat", y_id, file_dir=data_dir)
  do k = 1, size(lt)
    do i = 1, size(mt)
      write(y_id, trim(io_str)) y(:, i, k)
    end do
  end do
  close(y_id)

  call safe_open("n.dat", srn_id, file_dir=data_dir)
  do k = 1, size(lt)
    do i = 1, size(mt)
      write(srn_id, trim(io_str)) sr_n_0(:, i, k, 1)
    end do
  end do
  close(srn_id)

  call safe_open("sr1.dat", sr1_id, file_dir=data_dir)
  do k = 1, size(lt)
    do i = 1, size(mt)
      write(sr1_id, trim(io_str)) sr1(:, i, k)
    end do
  end do
  close(sr1_id)

  call safe_open("sr2.dat", sr2_id, file_dir=data_dir)
  do k = 1, size(lt)
    do i = 1, size(mt)
      write(sr2_id, trim(io_str)) sr2(:, i, k)
    end do
  end do
  close(sr2_id)

  call safe_open("sr3.dat", sr3_id, file_dir=data_dir)
  do k = 1, size(lt)
    do i = 1, size(mt)
      write(sr3_id, trim(io_str)) sr3(:, i, k)
    end do
  end do
  close(sr3_id)

end program find_solutions_initial
