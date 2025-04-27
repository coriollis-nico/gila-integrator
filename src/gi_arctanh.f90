! module gi_arctanh
!   !! TODO Write description
!   use iso_fortran_env,  &
!     only: sp => real32, &
!     dp => real64, &
!     qp => real128
!   use finite_diff
!   implicit none
!
!   private
!   public gi_arctanh_conditions
!   public gi_arctanh_friedmann, gi_arctanh_solution
!   public gi_arctanh_sr0, gi_arctanh_sr
!
!   type :: gi_arctanh_conditions
!     !! TODO
!     real(qp) :: Omega_M = 0.31_qp
!     !! \( \Omega_m \) initial matter density
!     real(qp) :: Omega_R = 8.4e-5_qp
!     !! \( \Omega_r \) initial radiation density
!     real(qp) :: Omega_dark = 0.69_qp
!     !! \( \Omega_{\Lambda} \) initial dark energy density
!     real(qp) :: l = 1.0e-27_qp
!     !! \( l = L H_0 \) early universe energy scale
!   end type
!
! contains
!
!   function gi_arctanh_friedmann(x, y, user_conditions)
!
!     !! TODO
!
!     real(qp), intent(in) :: x
!     !! \( x = \ln{\frac{a}{a_0}} \).
!     real(qp), intent(in) :: y
!     !! \( y = \frac{H}{H_0} \).
!     type(gi_arctanh_conditions), intent(in), optional :: user_conditions
!     !! User conditions for integration
!     real(qp) :: gi_arctanh_friedmann
!     ! Function output
!
!     type(gi_arctanh_conditions) :: cond
!     !! ODE parameters
!
!     if (present(user_conditions)) then
!       cond = user_conditions
!     end if
!
!     gi_arctanh_friedmann = - (1.0_qp/y) &
!                     * ( 1.0_qp - (cond%l * y)**4 ) &
!                     * ( (3.0_qp * cond%Omega_M / exp(3.0_qp * x) ) &
!                       + ( 4.0_qp * cond%Omega_R / exp(4.0_qp * x) ) )
!
!   end function gi_arctanh_friedmann
!
!   function gi_arctanh_solution(x, y0, user_conditions)
!
!     !! TODO
!
!     real(qp), intent(in) :: x(:)
!     !! [[gi_arctanh_friedmann:x]] range.
!     real(qp), intent(in) :: y0
!     !! initial [[gi_arctanh_friedmann:y]] value, corresponding to `x(1)`
!     type(gi_arctanh_conditions), intent(in), optional :: user_conditions
!
!     real(qp), dimension(size(x, 1)) :: gi_arctanh_solution
!
!     type(gi_arctanh_conditions) :: cond
!     real(qp) :: k1, k2, k3, k4
!     real(qp) :: h
!     integer :: i
!
!     if (present(user_conditions)) then
!       cond = user_conditions
!     end if
!
!     gi_arctanh_solution(1) = y0
!
!     do i = 2, size(x, 1)
!       h = x(i) - x(i-1)
!
!       k1 = gi_arctanh_friedmann(x(i-1), gi_arctanh_solution(i-1), cond)
!       k2 = gi_arctanh_friedmann(x(i-1) + h/2.0_qp, gi_arctanh_solution(i-1) + h*k1/2.0_qp, cond)
!       k3 = gi_arctanh_friedmann(x(i-1) + h/2.0_qp, gi_arctanh_solution(i-1) + h*k2/2.0_qp, cond)
!       k4 = gi_arctanh_friedmann(x(i-1) + h, gi_arctanh_solution(i-1) + h*k3, cond)
!
!       gi_arctanh_solution(i) = gi_arctanh_solution(i-1) &
!                               & + h * ( k1 + 2.0_qp * k2 + 2.0_qp * k3 + k4 ) / 6.0_qp
!     end do
!
!   end function gi_arctanh_solution
!
!   function gi_arctanh_sr0(x, y, l)
!
!     !! Returns a \( (N, \epsilon_0) \) two-column array from a [[gi_arctanh_solution]] curve.
!     !! Finds the `x` for which \( \bar{H} = l^{-1} \) to use as inflation exit point.
!
!     real(qp), intent(in) :: x(:)
!     !! \( \ln{\frac{a}{a_0}} \) equally spaced array.
!     real(qp), intent(in) :: y(:)
!     !! \( \frac{H}{H_0}} \) solution array.
!     real(qp), intent(in) :: l
!     !! Energy scale
!
!     real(qp) :: xi
!     !! `a_tilde` value for inflation exit
!     integer :: i_index
!     !! `a_tilde_i` index
!     real(qp) :: yi
!     !! `Hbar` value for `a_tilde_i`
!
!     integer :: i
!     ! iterator
!
!     real(qp), dimension(size(x, 1), 2) :: gi_arctanh_sr0
!     ! Output
!
!     if (size(x, 1) /= size(y, 1)) then
!       error stop "Arrays not of equal length"
!     end if
!
!     ! Finding a_tilde_i and Hbar_i
!     do i = 1, size(x, 1)
!       if ( y(i) >= 1._qp/l ) then
!         i_index = i
!         xi = x(i_index)
!         yi = y(i_index)
!         exit
!       end if
!     end do
!
!     do i = 1, size(x,1)
!       gi_arctanh_sr0(i,1) = x(i) - xi
!       gi_arctanh_sr0(i,2) = yi/y(i)
!     end do
!
!   end function gi_arctanh_sr0
!
!   function gi_arctanh_sr(x, y)
!
!     !! Given nth slow-roll data, finds \( \epsilon_{n+1} \). Assumes `\( x = \ln{a}{a_i} \)
!     !! is equally spaced.
!
!     real(qp), intent(in) :: x(:)
!     !! `x` output from [[gi_arctanh_sr]]
!     real(qp), intent(in) :: y(:)
!     !! \( \epsilon_n \)
!
!     real(qp), dimension(size(x, 1)) :: aux_array
!
!     integer :: i
!
!     real(qp), dimension(size(x, 1)) :: gi_arctanh_sr
!
!     if (size(x, 1) /= size(y, 1)) then
!       error stop "Arrays not of equal length"
!     end if
!
!     aux_array = [( log( abs( y(i) ) ), i = 1, size(x, 1) )]
!
!     gi_arctanh_sr = fd_c5curve(aux_array, x(2) - x(1), 1)
!
!   end function gi_arctanh_sr
!
! end module gi_arctanh
