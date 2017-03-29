module energy_eqns
  implicit none
  private
  public :: energy, temp
CONTAINS
  function energy(x,D,a) result(Eng)
    implicit none
    integer, parameter :: rk = 8
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x
    integer(ik) :: xsize, xs2
    real(rk) :: K, rho, b, ma, Eng
    real(rk), dimension(:) :: D, a
    real(rk), dimension(size(D)) :: xplus, xminus
    K = 0.025_rk
    rho = 2._rk
    b = 0.35_rk
    ma = 0.031_rk
    ! b = 0.35_rk
    ! K = 0.06_rk
    ! rho = 1._rk
    ! ma = 300._rk
    ! D = 0.03_rk
    ! a = 4.5_rk
    xsize = size(x)
    xs2 = xsize/2_ik
    xplus(2:xs2) = x(2:xs2) + x(1:xs2-1)
    xplus(1) = x(1) + x(xs2)

    xminus(2:xs2) = x(2:xs2) - x(1:xs2-1)
    xminus(1) = x(1) - x(xs2)
    Eng = (1._rk/(2._rk*ma)) * norm2(x(xs2+1:xsize))**2 + &
            sum(D*(exp(-a*x(1:xs2))-1)**2) + &
            sum((K/2._rk)*(1._rk + rho * exp(-b*xplus))*(xminus**2))
  end function energy
  function temp(x) result(T)
    implicit none
    integer, parameter :: rk = 8
    integer, parameter :: ik = 8
    real(rk), dimension(:), intent(in) :: x
    real(rk) :: T, kb=0.00008617_rk, ma=0.031_rk
    T = (sum(x(size(x)/2 + 1:size(x))**2))/((size(x)/2._rk) * kb * ma)
  end function temp
end module energy_eqns
