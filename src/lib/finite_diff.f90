module finite_diff
!! Module for calculating 5-point finite difference derivatives in
!! equally spaced curves
  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128
  use ieee_arithmetic,  &
    only: ieee_Value, ieee_quiet_nan
  implicit none

  private
  public fd_c5curve

contains

function center5_diff(y, h, n)

  !! Calculates \( y^{(n)} (x_0) \) in the middle point of dimension 5 array `y`

  real(qp), intent(in), dimension(5) :: y(-2:2)
  !! \( y = [y_{-2}. y_{-1}, y_0, y_1, y_2] \)
  real(qp), intent(in) :: h
  !! \( h = \Delta x \) (\(x_i\) equally spaced)
  integer, intent(in) :: n
  !! order of differentiation (1 or 2)
  real(qp) :: center5_diff

  select case(n)
    case(1)
      center5_diff = ( -y(2) + 8.0_qp * y(1) - 8.0_qp * y(-1) + y(-2) ) &
                   & / ( 12.0_qp * h )
    case(2)
      center5_diff = ( -y(2) + 16.0_qp * y(1) - 30.0_qp * y(0) &
                  & + 16.0_qp * y(-1) - y(-2) ) / ( 12.0_qp * h * h )
     case(3)
       center5_diff = ( y(2) -  2.0_qp * y(1) + 2.0_qp * y(-1) - y(-2) ) &
                   & / ( 2.0_qp * h**3 )
    case(4)
      center5_diff = ( y(2) - 4._qp*y(1) + 6._qp*y(0) - 4._qp*y(-1) + y(-2) ) &
                 & / (h**4)
    case default
      error stop "invalid n"
  end select

end function center5_diff

function fd_c5curve(y, h, n)

  !! Returns \( y^{(n)} (x) \)

  real(qp), dimension(:), intent(in) :: y
  real(qp), intent(in) :: h
  integer, intent(in) :: n

  real(qp), dimension( size(y) ) :: fd_c5curve

  integer :: i

  fd_c5curve(1:2) = ieee_value(fd_c5curve(1), ieee_quiet_nan)
  fd_c5curve(size(y)-1:size(y)) = ieee_value(fd_c5curve(1), ieee_quiet_nan)

  fd_c5curve(3:size(y)-2) = [(center5_diff(y(i-2:i+2), h, n), i = 3, size(y)-2)]

end function fd_c5curve

end module finite_diff
