module gi
  !! Provides interface to the GI Friedmann equation and an RK4 integration function
  !! (quadruple precision only)
  use iso_fortran_env,  &
    only: qp => real128
  use fdiff
  implicit none

  private
  public gi_conditions
  public write_gi_genconditions
  public gi_friedmann, gi_solution
  public gi_sr0, gi_sr

  type :: gi_conditions
    !! Derived type for specifying [[gi:gi_friedmann]] parameters
    real(qp) :: ld
    !! Early universe coefficient
    real(qp) :: l
    !! \( l = L H_0 \) early universe energy scale
    integer :: m
    !! TODO specify
    integer :: p
    !! TODO specify
    real(qp) :: Omega_M
    !! \( \Omega_m \) initial matter density
    real(qp) :: Omega_R
    !! \( \Omega_r \) initial radiation density
    real(qp) :: Omega_dark
    !! \( \Omega_{\Lambda} \) initial dark energy density
  end type

contains

  subroutine write_gi_genconditions(conditions, file_id)

    !! Writes relevant [[gi_conditions]] parameters to a file

    type(gi_conditions), intent(in) :: conditions
    integer, intent(in) :: file_id

    character(len=2) :: c = "# "
    !! Comment marker

    write(file_id, *) c//"Ω_M0 = ", conditions%Omega_M
    write(file_id, *) c//"Ω_R0 = ", conditions%Omega_R
    write(file_id, *) c//"Ω_Λ = ", conditions%Omega_dark
    write(file_id, *) c//"λ = ", conditions%ld
    write(file_id, *) c//"l = ", conditions%l
    write(file_id, *) c//"m = ", conditions%m
    write(file_id, *) c//"p = ", conditions%p


  end subroutine write_gi_genconditions

  pure function gi_friedmann(x, y, cond)

    !! Finds the value of \( \frac{d \bar{H}}{d \tilde{a}} \) of the GI Friedmann equations.

    real(qp), intent(in) :: x
    !! \( x = \ln{\frac{a}{a_0}} \).
    real(qp), intent(in) :: y
    !! \( y = \frac{H}{H_0} \).
    type(gi_conditions), intent(in) :: cond
    !! User conditions for integration
    real(qp) :: gi_friedmann
    ! Function output

    gi_friedmann = - ( 0.5_qp / ( y ) ) &
                     * ( 3.0_qp*cond%Omega_M/(exp(3.0_qp * x)) &
                        + 4.0_qp*cond%Omega_R/(exp(4.0_qp * x)) ) &
                     * ( ( 1._qp + cond%ld * cond%l**(2*cond%m - 2) &
                        * exp(cond%ld * cond%l**(2*cond%p)) ) &
                        / (1._qp &
                          + cond%ld * (cond%l * y)**(2*cond%m - 2) &
                            * (cond%m &
                                + cond%l * cond%p * (cond%l * y)**(2*cond%p)) &
                              * exp( cond%ld * (cond%l * y)**(2*cond%p) ) ) )

  end function gi_friedmann

  pure function gi_solution(x, y0, cond)

    !! Returns the RK4 solution to the [[gi_friedmann]] differential equation.

    ! TODO: Consider changing to 7 point deravative to reduce error

    real(qp), intent(in) :: x(:)
    !! [[gi_friedmann:x]] range.
    real(qp), intent(in) :: y0
    !! initial [[gi_friedmann:y]] value, corresponding to `x(1)`
    type(gi_conditions), intent(in) :: cond

    real(qp), dimension(size(x, 1)) :: gi_solution

    real(qp) :: k1, k2, k3, k4
    real(qp) :: h
    integer :: i

    gi_solution(1) = y0

    do i = 2, size(x, 1)
      h = x(i) - x(i-1)

      k1 = gi_friedmann(x(i-1), gi_solution(i-1), cond)
      k2 = gi_friedmann(x(i-1) + h/2.0_qp, gi_solution(i-1) + h*k1/2.0_qp, cond)
      k3 = gi_friedmann(x(i-1) + h/2.0_qp, gi_solution(i-1) + h*k2/2.0_qp, cond)
      k4 = gi_friedmann(x(i-1) + h, gi_solution(i-1) + h*k3, cond)

      gi_solution(i) = gi_solution(i-1) &
                    & + h * ( k1 + 2.0_qp * k2 + 2.0_qp * k3 + k4 ) / 6.0_qp
    end do

  end function gi_solution

  pure function gi_sr0(x, y, cond)

    !! Returns a \( (N, \epsilon_0) \) two-column array from a [[gi_solution]] curve.
    !! Finds the `x` for which \( \bar{H} = l^{-1} \) to use as inflation exit point.

    real(qp), intent(in) :: x(:)
    !! \( \ln{\frac{a}{a_0}} \) equally spaced array.
    real(qp), intent(in) :: y(:)
    !! \( \frac{H}{H_0}} \) solution array.
    type(gi_conditions), intent(in) :: cond

    real(qp) :: xi
    !! `a_tilde` value for inflation exit
    integer :: i_index
    !! `a_tilde_i` index
    real(qp) :: yi
    !! `Hbar` value for `a_tilde_i`

    integer :: i
    ! iterator

    real(qp), dimension(size(x, 1), 2) :: gi_sr0
    ! Output

    ! Finding a_tilde_i and Hbar_i
    do i = 1, size(x, 1)
      if ( y(i) >= 1._qp/cond%l ) then
        i_index = i
        xi = x(i_index)
        yi = y(i_index)
        exit
      end if
    end do

    do concurrent (i = 1:size(x,1))
      gi_sr0(i,1) = x(i) - xi
      gi_sr0(i,2) = yi/y(i)
    end do

  end function gi_sr0

  pure function gi_sr(x, y)

    !! Given nth slow-roll data, finds \( \epsilon_{n+1} \). Assumes `\( x = \ln{\frac{a}{a_i}} \)
    !! is equally spaced.

    real(qp), intent(in) :: x(:)
    !! `x` output from TODO
    real(qp), intent(in) :: y(:)
    !! \( \epsilon_n \)

    real(qp), dimension(size(x, 1)) :: gi_sr

    gi_sr = fdiff_c5curve(log(abs(y)), &
                              x(2) - x(1))

  end function gi_sr

end module gi
