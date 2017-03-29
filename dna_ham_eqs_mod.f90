module ham_eqs_mod
  implicit none
  private
  public :: ham_eqs_x, ham_eqs_p
CONTAINS
  function ham_eqs_x(x) result(dxdt)
    implicit none
    integer, parameter :: rk = 8
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x
    real(rk), dimension(size(x)/2) :: dxdt
    real(rk) :: ma
    integer(ik) :: xs, ii=0
    ! test compared to array comprehension
    ! ok, weird issue with array. might be faster, but bugs out
    ma = 0.031_rk
    ! ma = 300._rk
    xs = size(x)/2_ik
    do ii = 1,xs
      dxdt(ii) = x(xs+ii)/ma
    end do
  end function ham_eqs_x

  function ham_eqs_p(x,D,a) result(dpdt)
    ! takes vector x = (x1,x2,...,xN,p1,p2,...,pN)
    implicit none
    integer, parameter :: rk = 8
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x
    real(rk), dimension(size(x)/2) :: dpdt
    integer(ik) :: jj, xs2
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
    xs2 = size(x)/2_ik
    ! first (periodic bcs)
    dpdt(1) = D(1) * 2._rk * (-a(1) * exp(-a(1)*x(1))) * &
                    (exp(-a(1)*x(1)) - 1._rk) + &
                    (0.5_rk * K) * (-rho*b*exp(-b*(x(1)+x(xs2)))) * &
                    (x(1)-x(xs2))**2._rk + &
                    (0.5_rk * K) * (1._rk + rho * exp(-b*(x(1) + x(xs2)))) * &
                    2._rk * (x(1) - x(xs2)) + &
                    (0.5_rk * K) * (-rho*b*exp(-b*(x(2)+x(1)))) * &
                    (x(2)-x(1))**2 - &
                    (0.5_rk * K) * (1._rk + rho * exp(-b*(x(2) + x(1)))) * &
                    2._rk * (x(2) - x(1))


    ! middle section
    do jj = 2,xs2-1
      dpdt(jj) = D(jj) * 2._rk * (-a(jj) * exp(-a(jj)*x(jj))) * &
                      (exp(-a(jj)*x(jj)) - 1._rk) + &
                      (0.5_rk * K) * (-rho*b*exp(-b*(x(jj)+x(jj-1)))) * &
                      (x(jj)-x(jj-1))**2 + &
                      (0.5_rk * K) * (1._rk + rho * exp(-b*(x(jj) + x(jj-1)))) * &
                      2._rk * (x(jj) - x(jj-1)) + &
                      (0.5_rk * K) * (-rho*b*exp(-b*(x(jj+1)+x(jj)))) * &
                      (x(jj+1)-x(jj))**2 - &
                      (0.5_rk * K) * (1._rk + rho * exp(-b*(x(jj+1) + x(jj)))) * &
                      2._rk * (x(jj+1) - x(jj))
    end do

    ! last
    dpdt(xs2) = D(xs2) * 2._rk * (-a(xs2) * exp(-a(xs2)*x(xs2))) * &
                    (exp(-a(xs2)*x(xs2)) - 1._rk) + &
                    (0.5_rk * K) * (-rho*b*exp(-b*(x(xs2)+x(xs2-1)))) * &
                    (x(xs2)-x(xs2-1))**2 + &
                    (0.5_rk * K) * (1._rk + rho * exp(-b*(x(xs2) + x(xs2-1)))) * &
                    2._rk * (x(xs2) - x(xs2-1)) + &
                    (0.5_rk * K) * (-rho*b*exp(-b*(x(1)+x(xs2)))) * &
                    (x(1)-x(xs2))**2 - &
                    (0.5_rk * K) * (1._rk + rho * exp(-b*(x(1) + x(xs2)))) * &
                    2._rk * (x(1) - x(xs2))
  end function ham_eqs_p
end module ham_eqs_mod
