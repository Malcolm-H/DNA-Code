!!! not up to date


!!!
module ham_eqs_mod
  implicit none
  private
  public :: ham_eqs
CONTAINS
  function ham_eqs(x,atgc) result(dxdt)
    ! takes vector x = (x1,x2,...,xN,p1,p2,...,pN)
    implicit none
    integer, parameter :: rk = 16
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x, atgc ! atgc: 0 for at, 1 for gc
    real(rk), dimension(size(x)) :: dxdt
    integer(ik) :: ii, jj, kk, xsize, xs2
    real(rk) :: K=0.025_rk, rho=2._rk, b=0.35_rk
    real(rk), dimension(size(atgc),size(atgc)) :: Dgc, Dat, agc, aat
    real(rk), dimension(:) :: D, a
    xsize = size(x)
    xs2 = xsize/2_ik
    ! coefficient vectors depending on realisation
    do ii = 1,xsize/2_ik
      Dgc(ii,ii) = 0.075_rk
      Dat(ii,ii) = 0.05_rk
      agc(ii,ii) = 6.9_rk
      aat(ii,ii) = 4.2_rk
    end do
    D = MATMUL(Dgc,atgc) + MATMUL(Dat, 1._rk - atgc)
    a = MATMUL(agc,atgc) + MATMUL(aat, 1._rk - atgc)

    ! first (periodic bcs)
    dxdt(1) = x(xs2+1_ik)
    dxdt(xs2+1_ik) = D(1) * 2._rk * (-a(1) * exp(-a(1)*x(xsize))) * &
                    (exp(-a(1)*x(1)) - 1._rk) + &
                    (K/2._rk) * (-rho*b*exp(-b*(x(1)+x(xsize)))) * &
                    (x(1)-x(xsize))**2._rk + &
                    (K/2._rk) * (1._rk + rho * exp(x(1) + x(xsize))) * &
                    2 * (x(1) - x(xsize))
    ! middle N-2
    do jj = 2,xs2-1
      dxdt(jj) = x(jj + xs2)
      dxdt(jj+xs2) = D(jj) * 2._rk * (-a(jj) * exp(-a(jj)*x(jj))) * &
                      (exp(-a(jj)*x(jj)) - 1._rk) + &
                      (K/2._rk) * (-rho*b*exp(-b*(x(jj)+x(jj-1)))) * &
                      (x(jj)-x(jj-1))**2._rk + &
                      (K/2._rk) * (1._rk + rho * exp(x(jj) + x(jj-1))) * &
                      2 * (x(jj) - x(jj-1))
    end do

    ! last one
    dxdt(xs2) = x(xsize)
    dxdt(xsize) = D(xsize) * 2._rk * (-a(xsize) * exp(-a(xsize)*x(1))) * &
                    (exp(-a(xsize)*x(xsize)) - 1._rk) + &
                    (K/2._rk) * (-rho*b*exp(-b*(x(xsize)+x(1)))) * &
                    (x(xsize)-x(1))**2._rk + &
                    (K/2._rk) * (1._rk + rho * exp(x(xsize) + x(1))) * &
                    2 * (x(xsize) - x(1))
    end function ham_eqs
end module ham_eqs_mod
