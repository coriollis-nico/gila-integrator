program fd_test
   !! Test for [[finite_diff]] module
   use iso_fortran_env,  &
      only: qp => real128
   use finite_diff
   implicit none

   integer :: j

   integer, parameter :: n = 100
   real(qp), parameter :: h = 1.0_qp / n
   real(qp), dimension(n), parameter :: x = [( j*h, j = 1, n )]
   real(qp), dimension(n), parameter :: y = [( sin(x(j)), j = 1, n )]
   real(qp), dimension(n), parameter :: w1 = [( cos(x(j)), j = 1, n )]
   real(qp), dimension(n), parameter :: w2 = [( - sin(x(j)), j = 1, n )]
   real(qp), dimension(n), parameter :: w3 = [( - cos(x(j)), j = 1, n )]
   real(qp), dimension(n) :: z1
   real(qp), dimension(n) :: z2
   real(qp), dimension(n) :: z3
   real(qp), dimension(n) :: z4

! ------------------------------------------------------------------------------------ !

   z1 = fd_c5curve(y, h, 1)
   z2 = fd_c5curve(y, h, 2)
   z3 = fd_c5curve(y, h, 3)
   z4 = fd_c5curve(y, h, 4)

   print *, "Testing 1st derivative of sin(x)"
   do j = 3, n-2
      call diff_test( z1(j), w1(j) )
   end do
   print *, "Ok"

   print *, "Testing 2nd derivative of sin(x)"
   do j = 3, n-2
      call diff_test( z2(j), w2(j) )
   end do
   print *, "Ok"

   print *, "Testing 3rd derivative of sin(x)"
   do j = 3, n-2
      call diff_test( z3(j), w3(j) )
   end do
   print *, "Ok"

   print *, "Testing 4th derivative of sin(x)"
   do j = 3, n-2
      call diff_test( z4(j), y(j) )
   end do
   print *, "Ok"

contains

   function rel_error(a, b)

      real(qp), intent(in) :: a
      real(qp), intent(in) :: b

      real(qp) :: rel_error

      rel_error = ( a - b )/( b )

   end function rel_error

   subroutine diff_test(a, b)

      real(qp), intent(in) :: a
      real(qp), intent(in) :: b

      real(qp), parameter :: max_error = 0.001_qp

      if ( abs( rel_error(a,b) ) >= max_error ) then
         error stop "Difference is too big"
      end if

   end subroutine diff_test

end program fd_test
