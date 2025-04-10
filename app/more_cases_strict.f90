program more_cases
  !! For specified [[gila_conditions:m]], [[gila_conditions:p]], [[gila_conditions:p]]
  !! combinations calculates \( y \) and \( \epsilon_{1 \to 4} \) slow-roll parameters.
  !! Finds if sr curves touch regions of interest.
  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128,&
          stdout => output_unit
  use gila
  use lib_io
  implicit none

  integer :: i, k, j, l_c, n_c, e_c

  integer, dimension(100), parameter  :: mt = [(i, i = 3, 102)]
  !! [[gila_conditions:m]]
  integer, dimension(100), parameter  :: pt = [(i, i = 1, 100)]
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

  real(qp), parameter :: e1_max = 0.0097_qp
  real(qp), parameter :: e2_min = 0.024_qp, e2_max = 0.041_qp
  real(qp), parameter :: e3_min = -0.34_qp, e3_max = 0.74_qp

  logical, dimension(size(mt), size(pt), size(lt), 3) :: in_sr_range

! ------------------------------------------------------------------------------------ !

conditions%cosmos="early"
conditions%beta=0.0_qp
conditions%l_tilde=0.0_qp
conditions%r=0
conditions%s=0
conditions%Omega_M = matter_density
conditions%Omega_R = rad_density
conditions%Omega_dark = dark_density

in_sr_range(:,:,:,:) = .true.

do k = 1, size(lt)

  conditions%l=lt(k)

  do i = 1, size(mt)
  do j = 1, size(pt)

    conditions%m=mt(i)
    conditions%p=pt(j)

    y = gila_solution(x, yi, conditions)

    slowroll_data(:,1:2) = slowroll0(x, y, conditions%l)
    do l_c = 3, 5
      slowroll_data(:,l_c) = slowroll(slowroll_data(:,1), slowroll_data(:,l_c-1))
    end do


    do n_c = 1, size(slowroll_data, 1)
      if ((-60.0_qp <= slowroll_data(n_c, 1)) .and. (slowroll_data(n_c, 1) <= -50.0_qp)) then
        if (slowroll_data(n_c, 3) >= e1_max) then
          in_sr_range(i, j, k, 1) = .false.
          exit
        end if
      end if
    end do

    do n_c = 1, size(slowroll_data, 1)
      if ((-60.0_qp <= slowroll_data(n_c, 1)) .and. (slowroll_data(n_c, 1) <= -50.0_qp)) then
        if ((slowroll_data(n_c, 4) < e2_min) .or. (slowroll_data(n_c, 4) > e2_max)) then
          in_sr_range(i, j, k, 2) = .false.
          exit
        end if
      end if
    end do

    do n_c = 1, size(slowroll_data, 1)
      if ((-60.0_qp <= slowroll_data(n_c, 1)) .and. (slowroll_data(n_c, 1) <= -50.0_qp)) then
        if ((slowroll_data(n_c, 5) < e3_min) .or. (slowroll_data(n_c, 5) > e3_max)) then
          in_sr_range(i, j, k, 3) = .false.
          exit
        end if
      end if
    end do

    print *, "m =", mt(i)
    print *, "p =", pt(j)
    print *, "l =", lt(k)
    print *, "Result:", in_sr_range(i, j, k, :)
    print *, "-------------"

  end do
  end do

end do

end program more_cases
