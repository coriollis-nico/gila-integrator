program sr_variance_normalized
  !! Finds the variance \( \sigma^2 \) of the \( \epsilon_i \) values wrt the PLANCK reported
  !! values.
  use iso_fortran_env,  &
    only: qp => real128, &
    stdout => output_unit
  use gi
  use lib_io
  implicit none

  integer :: i, k, j, l_c, n_c, n_x

  integer, dimension(30), parameter  :: mt = [(i, i = 3, 32)]
  !! [[gila_conditions:m]]
  integer, dimension(10), parameter  :: pt = [(j, j = 1, 10)]
  !! [[gila_conditions:p]]
  real(qp), dimension(3), parameter :: lt = [1.e-17_qp, 1.e-22_qp, 1.e-27_qp]
  !! [[gila_conditions:l]]

  real(qp), parameter :: xi = 0.0_qp
  !! Initial [[gila_friedmann:x]] value
  real(qp), parameter :: xf = -100.0_qp
  !! Final [[gila_friedmann:y]] value
  real(qp), parameter :: yi = 1.0_qp
  !! Initial [[gila_friedmann:y]] value
  integer, parameter :: n = 20000
  !! Number of integration steps

  real(qp), dimension(n), parameter :: x = [( n_x * (xf- xi)/n, n_x = 0, n-1 )]
  !! \( x \) values array
  real(qp), dimension(n, size(mt), size(pt), size(lt)) :: y
  !! Solutions array
  real(qp), dimension(n, size(mt), size(pt), size(lt), -1:3) :: sr_data
  !! \( \epsilon_i \) array
  real(qp), dimension(size(mt), size(pt), size(lt), 3) :: sr_nv

  real(qp), parameter :: matter_density = 0.31_qp
  real(qp), parameter :: rad_density = 8.4e-5_qp
  real(qp), parameter :: dark_density = 0.69_qp

  real(qp), parameter :: e1 = 0, e1_max = 0.0097_qp
  real(qp), parameter :: e2 = 0.032_qp
  real(qp), parameter :: e2_max = e2 + 0.009_qp
  real(qp), parameter :: e3 = 0.19_qp
  real(qp), parameter :: e3_max = e3 + 0.55_qp
  real(qp), parameter :: n_min = -60.0_qp, n_max = -30.0_qp

  character(len=*), parameter :: data_dir = "data/tab/sr_variance_normalized"
  integer :: v_id, vn_id

! ------------------------------------------------------------------------------------ !

!   call safe_open("variance.dat", v_id, file_dir=data_dir)
!   call safe_open("variance_normalized.dat", vn_id, file_dir=data_dir)

  sr_nv(:,:,:,:) = 0

  print *, "Finding solutions"

  do concurrent (k = 1:size(lt))
    do concurrent (i = 1:size(mt))
      do concurrent (j = 1:size(pt))

        y(:, i, j, k) = gi_solution(x, yi, gi_conditions(&
                                              ld = 1.0_qp, &
                                              Omega_M = matter_density, &
                                              Omega_R = rad_density, &
                                              Omega_dark = dark_density, &
                                              l = lt(k), &
                                              p = pt(j), &
                                              m = mt(i) &
                                              )&
                                      )

        print '(a, i3)', "m = ", mt(i)
        print '(a, i3)', "p = ", pt(j)
        print '(a, es7.1e2)', "l = ", lt(k)
        print '(a)', "-----------"

      end do
    end do
  end do

!   print *, "Finding N, ϵ0"
!
!   do concurrent (k = 1:size(lt))
!     do concurrent (i = 1:size(mt))
!       do concurrent (j = 1:size(pt))
!
!         y(:, i, j, k) = gila_solution(x, yi, gi_conditions(&
!                                               ld = 1.0_qp, &
!                                               Omega_M = matter_density, &
!                                               Omega_R = rad_density, &
!                                               Omega_dark = dark_density, &
!                                               l = lt(k), &
!                                               p = pt(j), &
!                                               m = mt(i) &
!                                               )&
!                                       )
!
!         print '(a, i3)', "m = ", mt(i)
!         print '(a, i3)', "p = ", pt(j)
!         print '(a, es7.1e2)', "l = ", lt(k)
!         print '(a)', "-----------"
!
!       end do
!     end do
!   end do


!         slowroll_data(:,1:2) = slowroll0(x, y, conditions%l)
!         do l_c = 3, 5
!           slowroll_data(:,l_c) = slowroll(slowroll_data(:,1), slowroll_data(:,l_c-1))
!         end do
!
!         avg_num01 = 0
!         avg_num02 = 0
!         avg_num03 = 0
!         avg_denom = 0
!         do n_c = 1, size(slowroll_data, 1)
!           if ((n_min <= slowroll_data(n_c, 1)) .and. (slowroll_data(n_c, 1) <= n_max)) then
!             avg_num01 = avg_num01 + (slowroll_data(n_c, 3) - e1)**2
!             avg_num02 = avg_num02 + (slowroll_data(n_c, 4) - e2)**2
!             avg_num03 = avg_num03 + (slowroll_data(n_c, 5) - e3)**2
!             avg_denom = avg_denom + 1
!           end if
!         end do
!
!         sr_diff_avg(i, j, k, 1) = (avg_num01/avg_denom)
!         sr_diff_avg(i, j, k, 2) = (avg_num02/avg_denom)
!         sr_diff_avg(i, j, k, 3) = (avg_num03/avg_denom)
!
!
!   write(v_id, '(a)') "m    p    l    de1    de2    de3"
!   do i = 1, size(mt)
!     do j = 1, size(pt)
!       do k = 1, size(lt)
!         write(v_id, *) mt(i), pt(j), lt(k), &
!           sr_diff_avg(i, j, k, 1), sr_diff_avg(i, j, k, 2), sr_diff_avg(i, j, k, 3)
!       end do
!     end do
!   end do
!   close(v_id)
!
!   write(vn_id, '(a)') "m    p    l    de1   de2   de3"
!   do i = 1, size(mt)
!     do j = 1, size(pt)
!       do k = 1, size(lt)
!         write(vn_id, *) mt(i), pt(j), lt(k), &
!           sr_diff_avg(i, j, k, 1)/((e1_max - e1)**2), &
!           sr_diff_avg(i, j, k, 2)/((e2_max - e2)**2), &
!           sr_diff_avg(i, j, k, 3)/((e3_max - e3)**2)
!       end do
!     end do
!   end do
!   close(vn_id)

end program sr_variance_normalized
