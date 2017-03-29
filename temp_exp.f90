module energy_eqn
  implicit none
  private
  public :: energy
CONTAINS
  function energy(x,atgc) result(Eng)
    implicit none
    integer, parameter :: rk = 8
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x, atgc ! atgc: 0 for at, 1 for gc
    integer(ik) :: ii, jj, kk, xsize, xs2
    real(rk) :: K=0.025_rk, rho=2._rk, b=0.35_rk, ma = 0.031_rk, Eng
    real(rk), dimension(size(atgc),size(atgc)) :: Dgc, Dat, agc, aat
    real(rk), dimension(size(atgc)) :: D, a, xplus, xminus
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
    xplus(2:xs2) = x(2:xs2) + x(1:xs2-1)
    xplus(1) = x(1) + x(xs2)

    xminus(2:xs2) = x(2:xs2) - x(1:xs2-1)
    xminus(1) = x(1) - x(xs2)
    Eng = (1._rk/(2._rk*ma)) * norm2(x(xs2+1:xsize))**2 +&
            sum(D*2._rk*(exp(-a*x(1:xs2))-1)**2) + &
            sum((K/2._rk)*(1._rk + rho * exp(-b*xplus))*(xminus**2))
  end function energy
end module energy_eqn

module ham_eqs_mod
  implicit none
  private
  public :: ham_eqs
CONTAINS
  function ham_eqs(x,atgc) result(dxdt)
    ! takes vector x = (x1,x2,...,xN,p1,p2,...,pN)
    implicit none
    integer, parameter :: rk = 8
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x, atgc ! atgc: 0 for at, 1 for gc
    real(rk), dimension(size(x)) :: dxdt
    integer(ik) :: ii, jj, kk, xsize, xs2
    real(rk) :: K=0.025_rk, rho=2._rk, b=0.35_rk, ma = 0.031_rk
    real(rk), dimension(size(atgc),size(atgc)) :: Dgc, Dat, agc, aat
    real(rk), dimension(size(atgc)) :: D, a
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
    dxdt(1) = x(xs2+1_ik)/ma
    dxdt(xs2+1_ik) = D(1) * 2._rk * (-a(1) * exp(-a(1)*x(1))) * &
                    (exp(-a(1)*x(1)) - 1._rk) + &
                    (K/2._rk) * (-rho*b*exp(-b*(x(1)+x(xs2)))) * &
                    (x(1)-x(xs2))**2._rk + &
                    (K/2._rk) * (1._rk + rho * exp(-b*(x(1) + x(xs2)))) * &
                    2._rk * (x(1) - x(xs2)) + &
                    (K/2._rk) * (-rho*b*exp(-b*(x(2)+x(1)))) * &
                    (x(2)-x(1))**2 - &
                    (K/2._rk) * (1._rk + rho * exp(-b*(x(2) + x(1)))) * &
                    2._rk * (x(2) - x(1))
    ! middle section
    do jj = 2,xs2-1
      dxdt(jj) = x(jj + xs2)/ma
      dxdt(jj+xs2) = D(jj) * 2._rk * (-a(jj) * exp(-a(jj)*x(jj))) * &
                      (exp(-a(jj)*x(jj)) - 1._rk) + &
                      (K/2._rk) * (-rho*b*exp(-b*(x(jj)+x(jj-1)))) * &
                      (x(jj)-x(jj-1))**2 + &
                      (K/2._rk) * (1._rk + rho * exp(-b*(x(jj) + x(jj-1)))) * &
                      2._rk * (x(jj) - x(jj-1)) + &
                      (K/2._rk) * (-rho*b*exp(-b*(x(jj+1)+x(jj)))) * &
                      (x(jj+1)-x(jj))**2 - &
                      (K/2._rk) * (1._rk + rho * exp(-b*(x(jj+1) + x(jj)))) * &
                      2._rk * (x(jj+1) - x(jj))
    end do

    ! last
    dxdt(xs2) = x(xsize)/ma
    dxdt(xsize) = D(xs2) * 2._rk * (-a(xs2) * exp(-a(xs2)*x(xs2))) * &
                    (exp(-a(xs2)*x(xs2)) - 1._rk) + &
                    (K/2._rk) * (-rho*b*exp(-b*(x(xs2)+x(xs2-1)))) * &
                    (x(xs2)-x(xs2-1))**2._rk + &
                    (K/2._rk) * (1._rk + rho * exp(-b*(x(xs2) + x(xs2-1)))) * &
                    2._rk * (x(xs2) - x(xs2-1)) + &
                    (K/2._rk) * (-rho*b*exp(-b*(x(1)+x(xs2)))) * &
                    (x(1)-x(xs2))**2 - &
                    (K/2._rk) * (1._rk + rho * exp(-b*(x(1) + x(xs2)))) * &
                    2._rk * (x(1) - x(xs2))
  end function ham_eqs
end module ham_eqs_mod

program testy
  use energy_eqn, only : energy
  use ham_eqs_mod, only: ham_eqs
  implicit none
  real(8) :: x(20), y(10), eng, dxdt(20)
  integer(8) :: ii=0, tname=10
  character(len=8) :: x1, fmt
  fmt = '(I2.2)'
  write (x1,fmt) tname
  do ii=1,10
    x(ii) = 2.00 * 1.000
    x(ii+10) = 1.000
    y(ii) = 1.00
  end do
  dxdt = ham_eqs(x,y)
  open(unit=12, file='testing'//'_this_'//trim(x1)//'.txt')
  write(12,*) 101
  close(12)
  print*, dxdt
  ! eng = energy(x,y)
  ! print*, eng
end program testy
