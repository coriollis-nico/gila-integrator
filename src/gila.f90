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
public save_gila_genconditions, save_gila_limconditions
public gila_friedmann, gila_friedmann_limit, gila_solution, gila_solution_limit

type :: gila_conditions
  !! Derived type for specifying [[gila:gila_friedmann]] parameters
  character(5) :: cosmos = "gila"
  !! For selecting [[gila:gila_friedmann]] variation. Possible values are
  !!
  !! - `"gila"` (default) for the standard evaluation
  !! - `"early"` for computing \( \beta = 0 \) efficiently
  !! - `"late"` for computing \( \lambda = 0 \) efficiently
  !! - `"gr"` for computing \( \lambda = \beta = 0 \) (general relativity) efficiently
  logical :: dark = .false.
  !! For specifying \( \Lambda \) equal (`.false.`) or unequal (`.true.`) to `0`.
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
  real(qp) :: a0bar = 1.0_qp
  !! Initial \( \bar{a} = \frac{a}{a0} \) value
  real(qp) :: H0bar = 1.0_qp
  !! Initial \( \bar{H} = \frac{H}{H0} \) value
  real(qp) :: Omega_M = 0.999916_qp
  !! \( \Omega_m \) initial matter density
  real(qp) :: Omega_R = 8.4e-5_qp
  !! \( \Omega_r \) initial radiation density
  real(qp) :: Omega_dark = 0.0_qp
  !! \( \Omega_{\Lambda} \) initial dark energy density
end type

type(gila_conditions), parameter :: gila_grdefault &
  & = gila_conditions(cosmos="gr", &
                      & dark=.true.
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

  write(file_id, *) c//"* Cosmos: "//trim(conditions%cosmos)
  write(file_id, *) c//"dark: ", conditions%dark
  write(file_id, *) c//"a0bar = ", conditions%a0bar
  write(file_id, *) c//"H0bar = ", conditions%H0bar
  write(file_id, *) c//"Ω_M0 = ", conditions%Omega_M
  write(file_id, *) c//"Ω_R0 = ", conditions%Omega_R
  if ( conditions%dark == .true. ) then
    write(file_id, *) c//"Ω_Λ = ", conditions%Omega_dark
  end if
  if ( (trim(conditions%cosmos) == "gila") .or. (trim(conditions%cosmos) == "early") ) then
    write(file_id, *) c//"λ = ", conditions%lambda
    write(file_id, *) c//"l = ", conditions%l
    write(file_id, *) c//"m = ", conditions%m
    write(file_id, *) c//"p = ", conditions%p
  end if
  if ( (trim(conditions%cosmos) == "gila") .or. (trim(conditions%cosmos) == "late") ) then
    write(file_id, *) c//"β = ", conditions%lambda
    write(file_id, *) c//"l tilde = ", conditions%l_tilde
    write(file_id, *) c//"r = ", conditions%r
    write(file_id, *) c//"s = ", conditions%s
  end if


end subroutine save_gila_genconditions

subroutine save_gila_limconditions(conditions, file_id)

  !! [[save_gila_genconditions]] for the \( m,p \to \infty \) case.
  !! @warining
  !! \( r,s \to \infty \) not yet implemented in code
  !! @endwarning

  type(gila_conditions), intent(in) :: conditions
  integer, intent(in) :: file_id

  character(len=2) :: c = "# "
  !! Comment marker

  write(file_id, *) c//"* Cosmos: "//trim(conditions%cosmos)
  write(file_id, *) c//"dark: ", conditions%dark
  write(file_id, *) c//"a0bar = ", conditions%a0bar
  write(file_id, *) c//"H0bar = ", conditions%H0bar
  write(file_id, *) c//"Ω_M0 = ", conditions%Omega_M
  write(file_id, *) c//"Ω_R0 = ", conditions%Omega_R
  if ( conditions%dark == .true. ) then
    write(file_id, *) c//"Ω_Λ = ", conditions%Omega_dark
  end if
  if (trim(conditions%cosmos) == "early") then
    write(file_id, *) c//"λ = ", conditions%lambda
    write(file_id, *) c//"l = ", conditions%l
    write(file_id, *) c//"m -> oo"
    write(file_id, *) c//"p -> oo"
  else
    error stop "r,s limits not implemented"
!     write(file_id, *) c//"β = ", conditions%lambda
!     write(file_id, *) c//"l tilde = ", conditions%l_tilde
!     write(file_id, *) c//"r -> oo"
!     write(file_id, *) c//"s -> oo"
  end if

end subroutine save_gila_limconditions

function aux_num(lambda, l, m, p)

  !! Auxiliary function for use in [[gila_friedmann]].
  !! Calculates \( \lambda l^{2m-2} \exp{ \lambda l^{2p} } \) and
  !! \( - \beta \tilde{l}^{2r-2} \exp{ - \beta \tilde{l}^{2s} } \).

  real(qp), intent(in) :: lambda
  real(qp), intent(in) :: l
  integer, intent(in) :: m
  integer, intent(in) :: p

  real(qp) :: aux_num
  ! output

  aux_num = lambda * l**(2*m-2) * exp(lambda * l**(2p))

end function

function aux_denom(lambda, l, m, p, x)

  !! Denominator equivalent of [[aux_num]]
  !! for \( \lambda (xl)^{2m-2} ( m + \lambda p [xl]^{2p} ) \exp{ \lambda [xl]^{2p} } \)
  !! corresponding with \( -\beta \) etc.

  real(qp), intent(in) :: lambda
  real(qp), intent(in) :: l
  real(qp), intent(in) :: x
  integer, intent(in) :: m
  integer, intent(in) :: p

  real(qp) :: aux_denom
  ! output

  aux_denom = lambda * (x*l)**(2*m-2) * ( m + lambda * p * (l*x)**(2*p) ) &
            & * exp( lambda * (l*x)**(2p) )

end function aux_denom

function gila_friedmann(x, y, user_conditions)

  !! Finds the value of \( \frac{d \bar{H}}{d \bar{a}} \) of the GILA Friedmann equations.

  real(qp), intent(in) :: x
  !! \( x = \frac{a}{a_0} \).
  real(qp), intent(in) :: y
  !! \( y = \frac{H}{H_0} \).
  type(gila_conditions), intent(in), optional :: user_conditions
  !! User conditions for integration
  real(qp) :: gila_friedmann
  ! Function output

  type(gila_conditions) :: cond
  !! ODE parameters

  !- aux. values
  real(qp) :: aux_coeff
  !! Auxiliary for \( - \frac{1}{2 \bar{a} \bar{H}} \)
  real(qp) :: aux_densities
  !! Auxiliary for densities factor
  real(qp) :: aux_param
  !! Auxiliary for parameter fraction

  if ( present(user_conditions) ) then
    cond = user_conditions
  end if

  aux_coeff = - 0.5_qp / ( x * y )

  select case(cond%dark)
  case(.true.)
    aux_densities = 3.0_qp*cond%Omega_M/(x**3) + 4.0_qp*cond%Omega_R/(x**4) - 3.0_qp*cond%Omega_dark
  case(.false.)
    aux_densities = 3.0_qp*cond%Omega_M/(x**3) + 4.0_qp*cond%Omega_R/(x**4)
  case default
    error stop "Invalid 'dark' value"
  end select

  select case(trim(cond%cosmos))
  case("gr")
    aux_param = 1.0_qp
  case("early")
    aux_param = ( 1.0_qp + aux_num(cond%lambda, cond%l, cond%m, cond%p) ) &
              & / ( 1.0_qp + aux_denom(cond%lambda, cond%l, cond%m, cond%p, x) )
  case("late")
    aux_param = ( 1.0_qp + aux_num(-cond%beta, cond%l_tilde, cond%r, cond%s) ) &
              & / ( 1.0_qp + aux_denom(-cond%beta, cond%l_tilde, cond%r, cond%s, x) )
  case("gila")
    aux_param = ( 1.0_qp + aux_num(cond%lambda, cond%l,       cond%m, cond%p) &
                         + aux_num(-cond%beta,  cond%l_tilde, cond%r, cond%s) ) &
              & / ( 1.0_qp + aux_denom(cond%lambda, cond%l,       cond%m, cond%p, x) &
                         & + aux_denom(-cond%beta,  cond%l_tilde, cond%r, cond%s, x) )
                          !! NOTE the - sign before cond%beta
    case default
      error stop "Invalid cosmos selection"
  end select

  gila_friedmann = aux_coeff * aux_densities * aux_param

end function gila_friedmann

function gila_friedmann_limit(x, y, user_conditions)

  !! Finds the value \( y'(x) \) of [[gila_friedmann]] when \( m,p \to \infty \)
  !!
  !! TODO: incorporate `beta != 0` case

  real(qp), intent(in) :: x
  !! [[gila_friedmann:x]]
  real(qp), intent(in) :: y
  !! [[gila_friedmann:y]]
  type(gila_conditions), intent(in), optional :: user_conditions
  !! [[gila_friedmann:user_conditions]]
  real(qp) :: gila_friedmann_limit
  ! Function output

  type(gila_conditions) :: cond
  !! Conditions used by function

  !- aux. values
  real(qp) :: aux_coeff
  !! Auxiliary for \( - \frac{1}{2 \bar{a} \bar{H}} \)
  real(qp) :: aux_densities
  !! Auxiliary for densities factor

  if ( present(user_conditions) ) then
    cond = user_conditions
  end if

  aux_coeff = - 0.5_qp / ( x * y )

  select case(cond%dark)
  case(.true.)
    aux_densities = 3.0_qp*cond%Omega_M/(x**3) + 4.0_qp*cond%Omega_R/(x**4) - 3.0_qp*cond%Omega_dark
  case(.false.)
    aux_densities = 3.0_qp*cond%Omega_M/(x**3) + 4.0_qp*cond%Omega_R/(x**4)
  case default
    error stop "Invalid 'dark' value"
  end select

  if ( cond%l >= 1.0_qp ) then
    error stop "l >= 1 not implemented"
  end if

  select case(trim(cond%cosmos))
  case("early")
    if ( y >= 1.0_qp/l ) then
      gila_friedmann_limit = 0.0_qp
    else if ( (0 < y) .and. (y < 1.0_qp/l) )
      gila_friedmann_limit = aux_coeff * aux_densities
    else
      error stop "Condition not implemented"
    end if
  case default
    error stop "Cosmos selection not implemented"
  end select

end function gila_friedmann_limit

function gila_solution(x, y0, n, user_conditions)

  !! Returns the RK4 solution to the [[gila:gila_friedmann]] differential equation, as
  !! well as up to \( \epsilon_n \). `x_i` MUST BE EQUALLY SPACED.
  !!
  !! For ´y = gila_solution(x, y0, user_conditions)´,
  !!
  !! - ´y(i,1)´ is \( y(x_i) \)
  !! - ´y(i,j)´ is \( \epsilon_j \)
  !! TODO: Consider changing to 7 point deravative to reduce error

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
  if (n >= 1) then
    gila_solution(1, 2) = gila_friedmann(x(1), y0, cond)
  end if

  do i = 2, size(x, 1)
    h = x(i) - x(i-1)

    k1 = gila_friedmann(x(i-1), gila_solution(i-1, 1), cond)
    k2 = gila_friedmann(x(i-1) + h/2.0_qp, gila_solution(i-1, 1) + h*k1/2.0_qp, cond)
    k3 = gila_friedmann(x(i-1) + h/2.0_qp, gila_solution(i-1, 1) + h*k2/2.0_qp, cond)
    k4 = gila_friedmann(x(i-1) + h, gila_solution(i-1, 1) + h*k3, cond)

    gila_solution(i, 1) = gila_solution(i-1, 1) &
                        & + h * ( k1 + 2.0_qp * k2 + 2.0_qp * k3 + k4 ) / 6.0_qp
    if (n >= 1) then
      gila_solution(i, 2) = -k1
    end if
  end do

  do i = 2, n
    gila_solution(:,i+1) = fd_c5curve(gila_solution(:,i), x(2)-x(1), 1)
    do j = 1, size(x,1)
      gila_solution(j,i+1) = gila_solution(j,i+1) / gila_solution(j,i)
    end do
  end do

end function gila_solution

function gila_solution_limit(x, y0, n, user_conditions)

  !! Returns the RK4 solution to the [[gila:gila_friedmann_limit]] differential equation, as
  !! well as up to \( \epsilon_n \).
  !!
  !! For ´y = gila_solution_limit(x, y0, user_conditions)´,
  !!
  !! - ´y(i,1)´ is \( y(x_i) \)
  !! - ´y(i,j)´ is \( \epsilon_j \)

  ! BUG: this is copy-pasted from gila_soultion. Surely a better solution exists...

  real(qp), intent(in) :: x(:)
  real(qp), intent(in) :: y0
  type(gila_conditions), intent(in), optional :: user_conditions
  integer, intent(in) :: n
  !! If `n == 0`, only finds the solution curve

  real(qp), dimension(size(x, 1), n+1 ) :: gila_solution_limit

  type(gila_conditions) :: cond
  real(qp) :: k1, k2, k3, k4
  real(qp) :: h
  integer :: i, j

  if (present(user_conditions)) then
    cond = user_conditions
  end if

  gila_solution_limit(1, 1) = y0
  gila_solution_limit(1, 2) = gila_friedmann_limit(x(1), y0, cond)

  do i = 2, size(x, 1)
    h = x(i) - x(i-1)

    k1 = gila_friedmann_limit(x(i-1), gila_solution_limit(i-1, 1), cond)
    k2 = gila_friedmann_limit(x(i-1) + h/2.0_qp, gila_solution_limit(i-1, 1) + h*k1/2.0_qp, cond)
    k3 = gila_friedmann_limit(x(i-1) + h/2.0_qp, gila_solution_limit(i-1, 1) + h*k2/2.0_qp, cond)
    k4 = gila_friedmann_limit(x(i-1) + h, gila_solution_limit(i-1, 1) + h*k3, cond)

    gila_solution_limit(i, 1) = gila_solution_limit(i-1, 1) &
                        & + h * ( k1 + 2.0_qp * k2 + 2.0_qp * k3 + k4 ) / 6.0_qp
    if (n >= 1) then
      gila_solution_limit(i, 2) = -k1
    end if
  end do

  do i = 2, n
    gila_solution_limit(:,i+1) = fd_c5curve(gila_solution_limit(:,i), x(2)-x(1), 1)
    do j = 1, size(x,1)
      gila_solution_limit(j,i+1) = gila_solution_limit(j,i+1) / gila_solution_limit(j,i)
    end do
  end do

end function gila_solution_limit

end module gila
