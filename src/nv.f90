module nv
!! TODO write
  use iso_fortran_env,  &
    only: qp => real128
  implicit none

  private
  public nv_nvar

contains

  pure function nv_var(x, nmin, nmax, y, y_ref)

    !! TODO write

    real(qp), intent(in), dimension(:) :: x
    !! TODO
    real(qp), intent(in), dimension(:) :: y
    !! \( y(x) \)
    real(qp), intent(in) :: nmin, nmax
    !! TODO
    real(qp), intent(in) :: y_ref
    !! TODO

    real(qp) :: nv_var

    integer :: i
    real(qp) :: num
    integer :: denom

    num = 0
    denom = 0

    do i = 1, size(x)
      if (( x(i) <= nmax ) .and. ( x(i) >= nmin )) then
        num = num + (y(i) - y_ref)**2
        denom = denom + 1
      end if
    end do

    nv_var = num/( 1._qp * denom )

  end function nv_var

  pure function nv_nvar(x, nmin, nmax, y, y_ref, y_err)

    !! TODO write

    real(qp), intent(in), dimension(:) :: x
    !! TODO
    real(qp), intent(in), dimension(:) :: y
    !! \( y(x) \)
    real(qp), intent(in) :: nmin, nmax
    !! TODO
    real(qp), intent(in) :: y_ref, y_err
    !! TODO

    real(qp) :: nv_nvar

    nv_nvar = nv_var(x, nmin, nmax, y, y_ref)/(y_err**2)

  end function nv_nvar

end module nv
