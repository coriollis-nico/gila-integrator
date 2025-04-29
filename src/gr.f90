module gr
  !! Provides interface to the GR Friedmann equation and an RK4 integration function
  !! (quadruple precision only)
  use iso_fortran_env,  &
    only: qp => real128
  use fdiff
  implicit none

  private
  public gr_conditions, gr_conditions_default
  public write_gr_genconditions
  public gr_friedmann, gr_solution
  public gr_sr0, gr_sr

  type :: gr_conditions
    !! Derived type for specifying [[gi:gr_friedmann]] parameters
    real(qp) :: Omega_M
    !! \( \Omega_m \) initial matter density
    real(qp) :: Omega_R
    !! \( \Omega_r \) initial radiation density
    real(qp) :: Omega_dark
    !! \( \Omega_{\Lambda} \) initial dark energy density
  end type

  type(gr_conditions), parameter :: gr_conditions_default &
                                    = gr_conditions(Omega_M = 0.3153_qp, &
                                                    Omega_R = 9.02e-5_qp, &
                                                    Omega_dark = 0.6847_qp)
  !! Default parameters for GR integration, defined por convenience

contains

  subroutine write_gr_genconditions(conditions, file_id)

    !! Writes relevant [[gr_conditions]] parameters to a file

    type(gr_conditions), intent(in) :: conditions
    integer, intent(in) :: file_id

    character(len=2) :: c = "# "
    !! Comment marker

    write(file_id, *) c//"Ω_M0 = ", conditions%Omega_M
    write(file_id, *) c//"Ω_R0 = ", conditions%Omega_R
    write(file_id, *) c//"Ω_Λ = ", conditions%Omega_dark

  end subroutine write_gr_genconditions

  pure function gr_friedmann(x, y, cond)

    !! Finds the value of \( \frac{d \bar{H}}{d \tilde{a}} \) of the GI Friedmann equations.

    real(qp), intent(in) :: x
    !! \( x = \ln{\frac{a}{a_0}} \).
    real(qp), intent(in) :: y
    !! \( y = \frac{H}{H_0} \).
    type(gr_conditions), intent(in) :: cond
    !! User conditions for integration
    real(qp) :: gr_friedmann
    ! Function output

    gr_friedmann = - ( 0.5_qp / ( y ) ) &
                     * ( 3.0_qp*cond%Omega_M/(exp(3.0_qp * x)) &
                        + 4.0_qp*cond%Omega_R/(exp(4.0_qp * x)) )

  end function gr_friedmann

  pure function gr_solution(x, y0, cond)

    !! Returns the RK4 solution to the [[gr_friedmann]] differential equation.

    ! TODO: Consider changing to 7 point deravative to reduce error

    real(qp), intent(in) :: x(:)
    !! [[gr_friedmann:x]] range.
    real(qp), intent(in) :: y0
    !! initial [[gr_friedmann:y]] value, corresponding to `x(1)`
    type(gr_conditions), intent(in) :: cond

    real(qp), dimension(size(x, 1)) :: gr_solution

    real(qp) :: k1, k2, k3, k4
    real(qp) :: h
    integer :: i

    gr_solution(1) = y0

    do i = 2, size(x, 1)
      h = x(i) - x(i-1)

      k1 = gr_friedmann(x(i-1), gr_solution(i-1), cond)
      k2 = gr_friedmann(x(i-1) + h/2.0_qp, gr_solution(i-1) + h*k1/2.0_qp, cond)
      k3 = gr_friedmann(x(i-1) + h/2.0_qp, gr_solution(i-1) + h*k2/2.0_qp, cond)
      k4 = gr_friedmann(x(i-1) + h, gr_solution(i-1) + h*k3, cond)

      gr_solution(i) = gr_solution(i-1) &
                    & + h * ( k1 + 2.0_qp * k2 + 2.0_qp * k3 + k4 ) / 6.0_qp
    end do

  end function gr_solution

  pure function gr_sr0(x, y, l)

    !! Returns a \( (N, \epsilon_0) \) two-column array from a [[gr_solution]] curve.
    !! Finds the `x` for which \( \bar{H} = l^{-1} \) to use as inflation exit point.

    real(qp), intent(in) :: x(:)
    !! \( \ln{\frac{a}{a_0}} \) equally spaced array.
    real(qp), intent(in) :: y(:)
    !! \( \frac{H}{H_0}} \) solution array.
    real(qp), intent(in) :: l

    real(qp) :: xi
    !! `a_tilde` value for inflation exit
    integer :: i_index
    !! `a_tilde_i` index
    real(qp) :: yi
    !! `Hbar` value for `a_tilde_i`

    integer :: i
    ! iterator

    real(qp), dimension(size(x, 1), 2) :: gr_sr0
    ! Output

    ! Finding a_tilde_i and Hbar_i
    do i = 1, size(x, 1)
      if ( y(i) >= 1._qp/l ) then
        i_index = i
        xi = x(i_index)
        yi = y(i_index)
        exit
      end if
    end do

    do concurrent (i = 1:size(x,1))
      gr_sr0(i,1) = x(i) - xi
      gr_sr0(i,2) = yi/y(i)
    end do

  end function gr_sr0

  pure function gr_sr(x, y)

    !! Given nth slow-roll data, finds \( \epsilon_{n+1} \). Assumes `\( x = \ln{\frac{a}{a_i}} \)
    !! is equally spaced.

    real(qp), intent(in) :: x(:)
    !! `x` output from TODO
    real(qp), intent(in) :: y(:)
    !! \( \epsilon_n \)

    real(qp), dimension(size(x, 1)) :: gr_sr

    gr_sr = fdiff_c5curve(log(abs(y)), &
                              x(2) - x(1))

  end function gr_sr

end module gr
