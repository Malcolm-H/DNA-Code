module var_eqs_mod
  implicit none
  private
  public :: var_eqs_x, var_eqs_p
CONTAINS
  function var_eqs_x(w) result(dxdt)
    implicit none
    integer, parameter :: rk = 8
    integer, parameter :: ik = 4
    real(rk), dimension(:), intent(in) :: w
    real(rk), dimension(size(w)/2) :: dxdt
    real(rk) :: ma
    integer(ik) :: xs, ii=0
    ma = 0.031_rk
    ! ma = 300._rk
    xs = size(w)/2_ik
    do ii=1,xs
      dxdt(ii) = w(ii+xs)/ma
    end do
  end function var_eqs_x

  function var_eqs_p(x,w,D,a) result(dpdt)
    ! takes vector w = (dx1d,dx2,...,dxN,dp1,dp2,...,dpN)
    implicit none
    integer, parameter :: rk = 8
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x,w
    real(rk), dimension(size(x)/2) :: dpdt
    integer(ik) :: ii, xs
    real(rk) :: K, rho, b, ma
    real(rk), dimension(:) :: D, a
    K = 0.025_rk
    rho = 2._rk
    b = 0.35_rk
    ma = 0.031_rk

    ! b = 0.35_rk
    ! K = 0.06_rk
    ! rho = 1._rk
    ! ma = 300._rk
    xs = size(x)/2_ik
    ! first
    dpdt(1) = (0.5_rk*K*(rho*b**2*exp(-b*(x(1)+x(xs))))*(x(1)-x(xs))**2 &
                -K*(1._rk+rho*exp(-b*(x(1)+x(xs))))) * w(xs) &
                +(2._rk*D(1)*a(1)**2*exp(-a(1)*x(1))*(2._rk*exp(-a(1)*x(1))-1._rk)&
                +0.5_rk*K*(rho*b**2*exp(-b*(x(1)+x(xs))))*(x(1)-x(xs))**2&
                +2._rk*K*(-rho*b*exp(-b*(x(1)+x(xs))))*(x(1)-x(xs))&
                +K*(1._rk + rho*exp(-b*(x(1)+x(xs))))&
                +0.5_rk*K*(rho*b**2*exp(-b*(x(2)+x(1))))*(x(2)-x(1))**2&
                +2._rk*K*rho*b*exp(-b*(x(2)+x(1)))*(x(2)-x(1))&
                +K*(1._rk+rho*exp(-b*(x(2)+x(1)))))*w(1)&
                +(rho*K*exp(-b*(x(2)+x(1)))*(0.5_rk*b**2*(x(2)-x(1))**2-1._rk)-K)*w(2)
    ! middle N-2
    do ii=2,xs-1
      dpdt(ii) = (0.5_rk*K*(rho*b**2*exp(-b*(x(ii)+x(ii-1))))*(x(ii)-x(ii-1))**2 &
                  -K*(1._rk+rho*exp(-b*(x(ii)+x(ii-1))))) * w(ii-1) &
                  +(2._rk*D(ii)*a(ii)**2*exp(-a(ii)*x(ii))*(2._rk*exp(-a(ii)*x(ii))-1._rk)&
                  +0.5_rk*K*(rho*b**2*exp(-b*(x(ii)+x(ii-1))))*(x(ii)-x(ii-1))**2&
                  +2._rk*K*(-rho*b*exp(-b*(x(ii)+x(ii-1))))*(x(ii)-x(ii-1))&
                  +K*(1._rk + rho*exp(-b*(x(ii)+x(ii-1)))) &
                  +0.5_rk*K*(rho*b**2*exp(-b*(x(ii+1)+x(ii))))*(x(ii+1)-x(ii))**2&
                  +2._rk*K*rho*b*exp(-b*(x(ii+1)+x(ii)))*(x(ii+1)-x(ii))&
                  +K*(1._rk+rho*exp(-b*(x(ii+1)+x(ii)))))*w(ii)&
                  +(rho*K*exp(-b*(x(ii+1)+x(ii)))*(0.5_rk*b**2*(x(ii+1)-x(ii))**2-1._rk)-K)*w(ii+1)
    end do
    ! last
    dpdt(xs) = (0.5_rk*K*(rho*b**2*exp(-b*(x(xs)+x(xs-1))))*(x(xs)-x(xs-1))**2 &
                -K*(1._rk+rho*exp(-b*(x(xs)+x(xs-1))))) * w(xs-1) &
                +(2._rk*D(xs)*a(xs)**2*exp(-a(xs)*x(xs))*(2._rk*exp(-a(xs)*x(xs))-1._rk)&
                +0.5_rk*K*(rho*b**2*exp(-b*(x(xs)+x(xs-1))))*(x(xs)-x(xs-1))**2&
                +2._rk*K*(-rho*b*exp(-b*(x(xs)+x(xs-1))))*(x(xs)-x(xs-1))&
                +K*(1._rk + rho*exp(-b*(x(xs)+x(1))))&
                +0.5_rk*K*(rho*b**2*exp(-b*(x(1)+x(xs))))*(x(1)-x(xs))**2&
                +2._rk*K*rho*b*exp(-b*(x(1)+x(xs)))*(x(1)-x(xs))&
                +K*(1._rk+rho*exp(-b*(x(1)+x(xs)))))*w(xs)&
                +(rho*K*exp(-b*(x(1)+x(xs)))*(0.5_rk*b**2*(x(1)-x(xs))**2-1._rk)-K)*w(1)
  end function var_eqs_p
end module var_eqs_mod
