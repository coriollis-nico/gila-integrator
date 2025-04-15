program sr_variance
   !! Finds the variance \( \sigma^2 \) of the \( \epsilon_i \) values wrt the PLANCK reported
   !! values.
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
   integer, dimension(100), parameter  :: pt = [(j, j = 1, 100)]
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
   real(qp), parameter :: e2 = 0.032_qp
   real(qp), parameter :: e2_min = e2 - 0.008_qp, e2_max = e2 + 0.009_qp
   real(qp), parameter :: e3 = 0.19_qp
   real(qp), parameter :: e3_min = e3 - 0.53_qp, e3_max = e3 + 0.55_qp

   real(qp) :: sr1_diff_avg(size(mt), size(pt), size(lt))
   real(qp) :: sr_diff_avg(size(mt), size(pt), size(lt), 2:3)

   real(qp) :: avg_num01, avg_num02, avg_num03, avg_denom

   character(len=*), parameter :: data_dir = "data/tab/sr_variance"
   integer :: out_id

! ------------------------------------------------------------------------------------ !

   call safe_open("variance.dat", out_id, file_dir=data_dir)

   conditions%cosmos="early"
   conditions%beta=0.0_qp
   conditions%l_tilde=0.0_qp
   conditions%r=0
   conditions%s=0
   conditions%Omega_M = matter_density
   conditions%Omega_R = rad_density
   conditions%Omega_dark = dark_density

   sr1_diff_avg(:, :, :) = 0
   sr_diff_avg(:,:,:,:) = 0

   do k = 1, size(lt)
      conditions%l=lt(k)
      do i = 1, size(mt)
         conditions%m=mt(i)
         do j = 1, size(pt)
            conditions%p=pt(j)

            print '(a, i3)', "m = ", mt(i)
            print '(a, i3)', "p = ", pt(j)
            print '(a, es7.1e2)', "l = ", lt(k)
            print '(a)', "-----------"

            y = gila_solution(x, yi, conditions)

            slowroll_data(:,1:2) = slowroll0(x, y, conditions%l)
            do l_c = 3, 5
               slowroll_data(:,l_c) = slowroll(slowroll_data(:,1), slowroll_data(:,l_c-1))
            end do

            avg_num01 = 0
            avg_num02 = 0
            avg_num03 = 0
            avg_denom = 0
            do n_c = 1, size(slowroll_data, 1)
               if ((-60.0_qp <= slowroll_data(n_c, 1)) .and. (slowroll_data(n_c, 1) <= -50.0_qp)) then
                  avg_num01 = avg_num01 + (slowroll_data(n_c, 3) - e1_max)
                  avg_num02 = avg_num02 + (slowroll_data(n_c, 4) - e2)**2
                  avg_num03 = avg_num03 + (slowroll_data(n_c, 5) - e3)**2
                  avg_denom = avg_denom + 1
               end if
            end do

            sr1_diff_avg(i, j, k) = avg_num01/avg_denom
            sr_diff_avg(i, j, k, 2) = avg_num02/avg_denom
            sr_diff_avg(i, j, k, 3) = avg_num03/avg_denom

         end do
      end do
   end do


   write(out_id, '(a)') "m    p    l    Δϵ1    σ²2    σ²3"
   do i = 1, size(mt)
      do j = 1, size(pt)
         do k = 1, size(lt)
            write(out_id, *) mt(i), pt(j), lt(k), &
               sr1_diff_avg(i, j, k), sr_diff_avg(i, j, k, 2), sr_diff_avg(i, j, k, 3)
         end do
      end do
   end do
   close(out_id)

end program sr_variance
