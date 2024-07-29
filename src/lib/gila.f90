module gila
  !! Provides interface to the GILA Friedmann equation and an RK4 integration function
  !! (quadruple precision only)
  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128
  use finite_diff
  implicit none

private
public gila_conditions, gila_grdefault
public save_gila_genconditions
public gila_friedmann, gila_solution

type :: gila_conditions
  !! Derived type for specifying [[gila:gila_friedmann]] parameters
  character(5):: cosmos = "gila"
  !! For selecting [[gila:gila_friedmann]] variation. Possible values are
  !!
  !! - `"gila"` (default) for the standard evaluation
  !! - `"early"` for computing \( \beta = 0 \) efficiently
  !! - `"late"` for computing \( \lambda = 0 \) efficiently
  !! - `"gr"` for computing \( \lambda = \beta = 0 \) (general relativity) efficiently
  real(qp) :: lambda = 1.0_qp
  !! Early universe coefficient
  real(qp) :: beta = 1.0_qp
  !! Late universe coefficient
  real(qp) :: l = 1.0e-27_qp
  !! \( l = L H_0 \) early universe energy scale
  real(qp) :: l_tilde = 1.0_qp
  !! \( \tilde{l} = \tilde{L} H_0 \) late universe energy scale
  integer :: m = 3
  !! TODO specify
  integer :: r = 3
  !! TODO specify
  integer :: p = 1
  !! TODO specify
  integer :: s = 1
  !! TODO specify
  real(qp) :: a0 = 1.0_qp
  !! Initial scale factor
  real(qp) :: omega_m = 0.999916_qp
  !! \( \Omega_m \) initial matter density
  real(qp) :: omega_r = 8.4e-5_qp
  !! \( \Omega_r \) initial radiation density
  real(qp) :: omega_de = 0.0_qp
  !! \( \Omega_{\Lambda} \) initial dark energy density
end type

type(gila_conditions), parameter :: gila_grdefault &
  & = gila_conditions(cosmos="gr", &
                      & lambda = 0.0_qp, beta = 0.0_qp,         &
                      & l = 0.0_qp, l_tilde = 0.0_qp,           &
                      & m = 0, p = 0, r = 0, s = 0,             &
                      & omega_m = 0.31_qp, omega_r = 8.4e-5_qp, &
                      & omega_de = 0.69_qp )
  !! Default parameters for GR integration, defined por convenience

contains

subroutine save_gila_genconditions(conditions, file_id)

  !! Saves relevant [[gila_conditions]] parameters to a file

  type(gila_conditions), intent(in) :: conditions
  integer, intent(in) :: file_id

  character(len=2) :: c = "# "
  !! Comment marker

  write(file_id, *) c//"- Cosmos: "//trim(conditions%cosmos)
  select case(conditions%cosmos)
    case("gr")
      write(file_id, *) c//"a0: ", conditions%a0
      write(file_id, *) c//"Omega_m: ", conditions%omega_m
      write(file_id, *) c//"Omega_r: ", conditions%omega_r
      write(file_id, *) c//"Omega_Lambda: ", conditions%omega_de
    case("early")
      write(file_id, *) c//"lambda: ", conditions%lambda
      write(file_id, *) c//"l: ", conditions%l
      write(file_id, *) c//"m: ", conditions%m
      write(file_id, *) c//"p: ", conditions%p
      write(file_id, *) c//"a0: ", conditions%a0
      write(file_id, *) c//"Omega_m: ", conditions%omega_m
      write(file_id, *) c//"Omega_r: ", conditions%omega_r
    case("late")
      write(file_id, *) c//"beta: ", conditions%beta
      write(file_id, *) c//"l_tilde: ", conditions%l_tilde
      write(file_id, *) c//"r: ", conditions%r
      write(file_id, *) c//"s: ", conditions%s
      write(file_id, *) c//"a0: ", conditions%a0
      write(file_id, *) c//"Omega_m: ", conditions%omega_m
      write(file_id, *) c//"Omega_r: ", conditions%omega_r
    case default
      write(file_id, *) c//"lambda: ", conditions%lambda
      write(file_id, *) c//"beta: ", conditions%beta
      write(file_id, *) c//"l: ", conditions%l
      write(file_id, *) c//"l_tilde: ", conditions%l_tilde
      write(file_id, *) c//"m: ", conditions%m
      write(file_id, *) c//"p: ", conditions%p
      write(file_id, *) c//"r: ", conditions%r
      write(file_id, *) c//"s: ", conditions%s
      write(file_id, *) c//"a0: ", conditions%a0
      write(file_id, *) c//"Omega_m: ", conditions%omega_m
      write(file_id, *) c//"Omega_r: ", conditions%omega_r
  end select

end subroutine save_gila_genconditions

function gila_friedmann(x, y, user_conditions)

  !! Finds the value \( y'(x) \) of the GILA Friedmann equations.
  !! Reduces to GR for \( \lambda = \beta = 0 \)

  real(qp), intent(in) :: x
  !! \( x = \ln{\left( \frac{a}{a_0} \right)} \).
  real(qp), intent(in) :: y
  !! \( y = \ln{\left( \frac{H}{H_0} \right)} \).
  type(gila_conditions), intent(in), optional :: user_conditions
  !! User conditions for integration
  real(qp) :: gila_friedmann
  !! Function output

  type(gila_conditions) :: cond
  !! Densities used by function

  if ( present(user_conditions) ) then
    cond = user_conditions
  end if

  select case(cond%cosmos)
    case("gr")
      gila_friedmann = -(0.5_qp) * ( 3.0_qp + cond%omega_r / exp(2.0_qp*y + 4.0_qp*x) &
      & - 3.0_qp * cond%omega_de )
    case("early")
      gila_friedmann = -(0.5_qp) * exp( -3.0_qp*x - 2.0_qp*y ) * cond%a0**(-3) &
      & * ( 3.0_qp * cond%omega_m + 4.0_qp * cond%omega_r / cond%a0 / exp(x) )             &
      & * ( 1.0_qp + cond%lambda * cond%l**(2 * cond%m - 2)                 &
      &     * exp( cond%lambda * cond%l**(2 * cond%p) ) )                   &
      & / ( 1.0_qp + cond%lambda * cond%l**(2 * cond%m - 2) * exp( (2 * cond%m - 2) * y ) &
      &       * exp( cond%lambda * cond%l**(2 * cond%p) * exp(2 * cond%p * y) )   &
      &       * ( cond%m + cond%lambda * cond%p * cond%l**(2 * cond%p)            &
      &       * exp( 2 * cond%p * y ) ) )
    case("late")
      gila_friedmann = -(0.5_qp) * exp( -3._qp*x - 2.0_qp*y ) * cond%a0**(-3) &
      & * ( 3.0_qp * cond%omega_m + 4.0_qp * cond%omega_r / cond%a0 / exp(x) )             &
      & * ( 1.0_qp - cond%beta   * cond%l_tilde**(2 * cond%r - 2)           &
      &     * exp( -cond%beta  * cond%l_tilde**(2 * cond%s) ) )             &
      & / ( 1.0_qp - cond%beta * cond%l_tilde**(2 * cond%r - 2) &
      &       * exp( (2 * cond%r - 2) * y ) &
      &       * exp( - cond%beta * cond%l_tilde**(2 * cond%s) &
      &       * exp(2 * cond%s * y) ) &
      &       * ( cond%r - cond%beta * cond%s * cond%l_tilde**(2 * cond%s) &
      &       * exp( 2 * cond%s * y ) ) )
    case("gila")
      gila_friedmann = -(0.5_qp) * exp( -3._qp*x - 2.0_qp*y ) * cond%a0**(-3) &
      & * ( 3.0_qp * cond%omega_m + 4.0_qp * cond%omega_r /cond%a0 / exp(x) )             &
      & * ( 1.0_qp + cond%lambda * cond%l**(2 * cond%m - 2)                 &
      &     * exp( cond%lambda * cond%l**(2 * cond%p) )                     &
      &            - cond%beta   * cond%l_tilde**(2 * cond%r - 2)           &
      &     * exp( -cond%beta  * cond%l_tilde**(2 * cond%s) ) )             &
      & / ( 1.0_qp + cond%lambda * cond%l**(2 * cond%m - 2) * exp( (2 * cond%m - 2) * y ) &
      &       * exp( cond%lambda * cond%l**(2 * cond%p) * exp(2 * cond%p * y) )   &
      &       * ( cond%m + cond%lambda * cond%p * cond%l**(2 * cond%p)            &
      &       * exp( 2 * cond%p * y ) ) &
      &            - cond%beta * cond%l_tilde**(2 * cond%r - 2) &
      &       * exp( (2 * cond%r - 2) * y ) &
      &       * exp( - cond%beta * cond%l_tilde**(2 * cond%s) &
      &       * exp(2 * cond%s * y) ) &
      &       * ( cond%r - cond%beta * cond%s * cond%l_tilde**(2 * cond%s) &
      &       * exp( 2 * cond%s * y ) ) )
    case default
      error stop "Invalid cosmos selection"
  end select

end function gila_friedmann

function gila_solution(x, y0, n, user_conditions)

  !! Returns the RK4 solution to the [[gila:gila_friedmann]] differential equation, as
  !! well as up to \( \epsilon_n \).
  !!
  !! For ´y = gila_solution(x, y0, user_conditions)´,
  !!
  !! - ´y(i,1)´ is \( y(x_i) \)
  !! - ´y(i,j)´ is \( \epsilon_j \)

  real(qp), intent(in) :: x(:)
  real(qp), intent(in) :: y0
  type(gila_conditions), intent(in), optional :: user_conditions
  integer, intent(in) :: n
  !! If `n == 0`, only finds the solution curve

  real(qp), dimension(size(x, 1), n+1 ) :: gila_solution

  type(gila_conditions) :: cond
  real(qp) :: k1, k2, k3, k4
  real(qp) :: h
  integer :: i, j

  if (present(user_conditions)) then
    cond = user_conditions
  end if

  gila_solution(1, 1) = y0
  gila_solution(1, 2) = gila_friedmann(x(1), y0, cond)

  do i = 2, size(x, 1)
    h = x(i) - x(i-1)

    k1 = gila_friedmann(x(i-1), gila_solution(i-1, 1), cond)
    k2 = gila_friedmann(x(i-1) + h/2.0_qp, gila_solution(i-1, 1) + h*k1/2.0_qp, cond)
    k3 = gila_friedmann(x(i-1) + h/2.0_qp, gila_solution(i-1, 1) + h*k2/2.0_qp, cond)
    k4 = gila_friedmann(x(i-1) + h, gila_solution(i-1, 1) + h*k3, cond)

    gila_solution(i, 1) = gila_solution(i-1, 1) &
                        & + h * ( k1 + 2.0_qp * k2 + 2.0_qp * k3 + k4 ) / 6.0_qp
    if (n >= 1) then
      gila_solution(i, 2) = k1
    end if
  end do

  do i = 2, n
    gila_solution(:,i+1) = fd_c5curve(gila_solution(:,2), x(2)-x(1), i-1)
    do j = 1, size(x,1)
      gila_solution(j,i+1) = gila_solution(j,i+1) / gila_solution(j,i)
    end do
  end do

end function gila_solution

end module gila
