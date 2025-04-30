program sr_normvar
  !! Finds the variance \( \sigma_n^2 \) of the \( \epsilon_i \) values wrt the PLANCK reported
  !! values.
  use iso_fortran_env,  &
    only: qp => real128, &
    stdout => output_unit
  use gi
  use io
  use nv
  implicit none

  integer :: i, k, j, e_c, n_x, this_div

  integer, parameter :: divs = 12
  integer, parameter :: msize = 20
  ! such that divs * msize = 240

  integer, dimension(divs*msize), parameter  :: mt = [(i, i = 3, divs*msize + 2)]
  !! [[gi_conditions:m]]
  integer, dimension(20), parameter  :: pt = [(j, j = 1, 20)]
  !! [[gi_conditions:p]]
  real(qp), dimension(3), parameter :: lt = [1.e-17_qp, 1.e-22_qp, 1.e-27_qp]
  !! [[gi_conditions:l]]

  real(qp), parameter :: xi = 0.0_qp
  !! Initial [[gi_friedmann:x]] value
  real(qp), parameter :: xf = -100.0_qp
  !! Final [[gi_friedmann:y]] value
  real(qp), parameter :: yi = 1.0_qp
  !! Initial [[gi_friedmann:y]] value
  integer, parameter :: n = 20000
  !! Number of integration steps

  real(qp), dimension(n), parameter :: x = [( n_x * (xf- xi)/n, n_x = 0, n-1 )]
  !! \( x \) values array
  real(qp), dimension(n, msize, size(pt), size(lt), -1:3) :: sr_data
  !! \( \epsilon_i \) array
  real(qp), dimension(msize, size(pt), size(lt), 3) :: sr_nv

  real(qp), parameter :: matter_density = 0.31_qp
  real(qp), parameter :: rad_density = 8.4e-5_qp
  real(qp), parameter :: dark_density = 0.69_qp

  real(qp), dimension(3), parameter :: sr_ref = [0._qp, 0.032_qp, 0.19_qp]
  real(qp), dimension(3), parameter :: sr_max = [0.0097_qp, &
                                                sr_ref(2) + 0.009_qp, &
                                                sr_ref(3) + 0.55_qp]
  real(qp), parameter :: n_min = -60.0_qp, &
                         n_max = -40.0_qp

  character(len=*), parameter :: data_dir = "data/tab/sr_normvar"
  integer :: vn_id

! ------------------------------------------------------------------------------------ !

  call safe_open("vn.dat", vn_id, file_dir=data_dir)
  write(vn_id, '(a)') "m    p    log_l    v1    v2    v3"

  do this_div = 0, divs-1

    print '(2(a, i2))', "Working ", this_div+1, "/", divs
    print '(a)', "! --------------------- !"
    print '(a)', "Finding N, ϵ0"

    do concurrent (k = 1:size(lt))
      do concurrent (j = 1:size(pt))
        do concurrent (i = 1 : msize)

          sr_data(:, i, j, k, -1:0) = gi_sr0(x, &
                                            gi_solution(x,&
                                              yi, &
                                              gi_conditions(&
                                                ld = 1.0_qp, &
                                                Omega_M = matter_density, &
                                                Omega_R = rad_density, &
                                                Omega_dark = dark_density, &
                                                l = lt(k), &
                                                p = pt(j), &
                                                m = mt(msize*this_div + i) &
                                              )&
                                              ), &
                                            gi_conditions(&
                                              ld = 1.0_qp, &
                                              Omega_M = matter_density, &
                                              Omega_R = rad_density, &
                                              Omega_dark = dark_density, &
                                              l = lt(k), &
                                              p = pt(j), &
                                              m = mt(msize*this_div + i) &
                                              ))

!           print '(a, i3)', "m = ", mt(msize*this_div + i)
!           print '(a, i3)', "p = ", pt(j)
!           print '(a, es8.1e2)', "l = ", lt(k)
!           print '(a)', "-----------"

        end do
      end do
    end do

    do e_c = 1, 3
      print '(a, i1)', "Finding ϵ_i:", e_c
      print '(a)', "--------------------"
      do concurrent (k = 1:size(lt))
        do concurrent (j = 1:size(pt))
          do concurrent (i = 1 : msize)

            sr_data(:, i, j, k, e_c) = gi_sr(x, sr_data(:, i, j, k, e_c - 1))

!             print '(a, i3)', "m = ", mt(msize*this_div + i)
!             print '(a, i3)', "p = ", pt(j)
!             print '(a, es8.1e2)', "l = ", lt(k)
!             print '(a)', "-----------"

          end do
        end do
      end do
    end do

    do e_c = 1, 3
      print '(a, i1)', "Finding variances:", e_c
      print '(a)', "--------------------"
      do concurrent (k = 1:size(lt))
        do concurrent (j = 1:size(pt))
          do concurrent (i = 1 : msize)

            sr_nv(i, j, k, e_c) = nv_nvar(sr_data(:, i, j, k, -1), &
                                  n_min, n_max, &
                                  sr_data(:, i, j, k, e_c), &
                                  sr_ref(e_c), sr_max(e_c))
!
!             print '(a, i3)', "m = ", mt(msize*this_div + i)
!             print '(a, i3)', "p = ", pt(j)
!             print '(a, es8.1e2)', "l = ", lt(k)
!             print '(a)', "-----------"

          end do
        end do
      end do
    end do

    do k = 1, size(lt)
    do j = 1, size(pt)
    do i = 1, msize
      write(vn_id, '(3(i5), 3(f30.25))') mt(msize*this_div + i), &
                                          pt(j), &
                                          int(log10(lt(k))), &
                                          sr_nv(i, j, k, :)
    end do
    end do
    end do

  end do

  close(vn_id)

end program sr_normvar
