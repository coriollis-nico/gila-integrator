module fdiff
!! Module for calculating 5-point finite difference derivatives in
!! equally spaced curves
  use iso_fortran_env,  &
    only: qp => real128
  implicit none

  private
  public fdiff_c5, fdiff_c5curve

contains

  pure function fdiff_c5(ym2, ym1, y1, y2, h)

    !! Calculates \( y'(x_0) \) in the middle point of dimension 5 array `y`

    real(qp), intent(in) :: ym2, ym1, y1, y2
    !! \( y = [y_{-2}. y_{-1}, y_0, y_1, y_2] \)
    real(qp), intent(in) :: h
    !! \( h = \Delta x \) (\(x_i\) equally spaced)
    real(qp) :: fdiff_c5

    fdiff_c5 = ( -y2 + 8.0_qp * y1 - 8.0_qp * ym1 + ym2 ) &
      & / ( 12.0_qp * h )

  end function fdiff_c5

  pure function fdiff_f3(y0, y1, y2, h)

    !! TODO write descr

    real(qp), intent(in) :: y0, y1, y2
    !! TODO
    real(qp), intent(in) :: h
    !! \( h = \Delta x \) (\(x_i\) equally spaced)
    real(qp) :: fdiff_f3

    fdiff_f3 = ( -3.0_qp*y0 + 4.0_qp*y1 - y2 ) &
      & / ( 2.0_qp * h )

  end function fdiff_f3

  pure function fdiff_b3(ym2, ym1, y0, h)

    !! TODO write descr

    real(qp), intent(in) :: ym2, ym1, y0
    !! TODO
    real(qp), intent(in) :: h
    !! \( h = \Delta x \) (\(x_i\) equally spaced)
    real(qp) :: fdiff_b3

    fdiff_b3 = ( 3.0_qp*y0 - 4.0_qp*ym1 + ym2 ) &
      & / ( 2.0_qp * h )

  end function fdiff_b3


  pure function fdiff_c5curve(y, h)

    !! Returns \( y'(x) \)

    real(qp), dimension(:), intent(in) :: y
    real(qp), intent(in) :: h

    real(qp), dimension( size(y) ) :: fdiff_c5curve

    integer :: i

    fdiff_c5curve(1) = fdiff_f3(y(1), y(2), y(3), h)
    fdiff_c5curve(2) = fdiff_f3(y(2), y(3), y(4), h)
    do i = 3, size(y)-2
      fdiff_c5curve(i) = fdiff_c5(y(i-2), y(i-1), y(i+1), y(i+2), h)
    end do
    fdiff_c5curve(size(y)-1) = fdiff_b3(y(size(y)-3), y(size(y)-2), y(size(y)-1), h)
    fdiff_c5curve(size(y)) = fdiff_b3(y(size(y)-2), y(size(y)-1), y(size(y)), h)

  end function fdiff_c5curve

end module fdiff
