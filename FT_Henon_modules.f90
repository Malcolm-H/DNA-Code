module hameq
  implicit none
  private
  public :: hameqs
CONTAINS
  function hameqs(x) result(dxdt)
    implicit none
    integer, parameter :: rk = 16
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x
    real(rk), dimension(size(x)) :: dxdt
    dxdt(1) = x(3)
    dxdt(2) = x(4)
    dxdt(3) = -(x(1) + 2._rk * x(1) * x(2))
    dxdt(4) = -(x(2) + x(1)**2._rk - x(2) ** 2._rk)
  end function hameqs
end module hameq
