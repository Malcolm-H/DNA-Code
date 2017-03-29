! most basic version of Ruth Forest order 4, not optimised

program Henon_Order4
  implicit none
  integer, parameter :: rk = 16
  integer, parameter :: ik = 8, ndof=4_ik, nfin=20_ik,nfoutx=19_ik,nfoute=18_ik
  integer, parameter :: nfoutl=17_ik, nfouts = 16_ik
  real(rk) :: h, E0, E, alphasum=0._rk, alpha1, alpha2, chi=0._rk, sali
  integer(ik) :: tf, niter, ii=1_ik, jj=1_ik, kk
  real(rk), dimension(ndof) :: x, veq, veq2, dxdt
  real(rk), dimension(4) :: c,d
  open(unit=nfin, file='FT_Henon.ini')
  read(nfin,*) E0
  read(nfin,*) x(1), x(2), x(4)
  read(nfin,*) veq
  read(nfin,*) veq2
  close(nfin)
  x(3) = sqrt(2._rk*E0 - x(1)**2._rk - x(2)**2._rk - x(4)**2._rk -2._rk * x(1)**2._rk * x(2) + (2.0_rk/3.0_rk) * x(2)**3_rk)
  h = 0.1_rk
  tf = 1000_ik
  niter = tf/h
  ! order 4 coefficients
  c(1) = 1._rk/(2._rk*(2._rk-2._rk**(1._rk/3._rk)))
  c(2) = (1._rk-2._rk**(1._rk/3._rk)) / (2._rk*(2._rk-2._rk**(1._rk/3._rk)))
  c(3) = c(2)
  c(4) = c(1)

  d(1) = 1._rk / (2._rk - 2._rk**(1._rk/3._rk))
  d(2) = - (2._rk**(1._rk/3._rk))/(2._rk - 2._rk**(1._rk/3._rk))
  d(3) = d(1)
  d(4) = 0._rk

  open(unit=nfoutx, file='Henon_Order4_xout.csv')
  open(unit=nfoute, file='Henon_Order4_Eout.csv')
  open(unit=nfoutl, file='Henon_Order4_lout.csv')
  open(unit=nfouts, file='Henon_Order4_sout.csv')

  write(nfoutx,*) x
  write(nfoute,*) 0._rk
  write(nfoutl,*) 1._rk/h
  write(nfouts,*) sqrt(2._rk)
  do ii=1_ik,niter
    do  jj=1_ik,size(c)
      dxdt = hameqs(x)
      ! x(1) = x(1) + h * c(jj) * x(3)
      ! x(2) = x(2) + h * c(jj) * x(4)
      x(1:2) = x(1:2) + h * c(jj) * dxdt(1:2)
      ! x(2) = x(2) + h * c(jj) * dxdt(2)
      veq(1) = veq(1) + h * c(jj) * veq(3)
      veq(2) = veq(2) + h * c(jj) * veq(4)
      veq2(1) = veq2(1) + h * c(jj) * veq2(3)
      veq2(2) = veq2(2) + h * c(jj) * veq2(4)
      dxdt = hameqs(x)
      ! x(3) = x(3) - h * d(jj) * (x(1) + 2._rk * x(1) * x(2))
      ! x(4) = x(4) - h * d(jj) * (x(2) + x(1)**2._rk - x(2) ** 2._rk)
      x(3) = x(3) + h * d(jj) * dxdt(3)
      x(4) = x(4) + h * d(jj)  * dxdt(4)
      veq(3) = veq(3) + d(jj) * h * ((-1._rk - 2._rk * x(2)) * veq(1) - 2._rk * x(1) * veq(2));
      veq(4) = veq(4) + d(jj) * h * (-2._rk * x(1) * veq(1) + (-1._rk + 2._rk * x(2)) * veq(2));
      veq2(3) = veq2(3) + d(jj) * h * ((-1._rk - 2._rk * x(2)) * veq2(1) - 2._rk * x(1) * veq2(2));
      veq2(4) = veq2(4) + d(jj) * h * (-2._rk * x(1) * veq2(1) + (-1._rk + 2._rk * x(2)) * veq2(2));
    end do
    alpha1 = norm2(veq)
    alpha2 = norm2(veq2)
    alphasum = alphasum + log(alpha1)
    chi = (1._rk/(ii*h)) * alphasum
    write(nfoutx,*) x
    call Eng(E,x)
    write(nfoute,*) abs((E-E0)/E0)
    write(nfoutl,*) chi
    do kk=1_ik,size(veq)
      veq(kk) = veq(kk)/alpha1
      veq2(kk) = veq2(kk)/alpha2
    end do
    sali = min(norm2(veq+veq2),norm2(veq-veq2))
    write(nfouts,*) sali
  end do
  close(nfoutx)
  close(nfoute)
  close(nfoutl)
  close(nfouts)
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

end program Henon_Order4

subroutine Eng(E,x)
  implicit none
  integer, parameter:: rk = 16
  real(rk) :: E
  real(rk), dimension(4) :: x
  E = 0.5_rk * (x(3)**2._rk + x(4)**2._rk) + 0.5_rk * (x(1)**2._rk + x(2)**2._rk) &
    + x(1) ** 2._rk * x(2) - (1._rk/3._rk) * x(2)**3._rk
end subroutine Eng
