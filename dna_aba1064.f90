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
    integer(ik) :: ii, jj, kk, xs2
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
    integer(ik) :: ii, jj, kk, xs
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
    integer(ik) :: ii, jj, kk, xsize, xs2
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

module rand_mod
  implicit none
  integer, parameter :: rk = 8
  integer, parameter :: ik = 8
CONTAINS
  ! from Francesco on StackOverflow
  subroutine init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
  end subroutine init_random_seed
end module rand_mod

program dnao4
  use ham_eqs_mod, only: ham_eqs_x, ham_eqs_p
  use energy_eqns, only: energy, temp
  use var_eqs_mod, only: var_eqs_x, var_eqs_p
  use rand_mod, only: init_random_seed
  implicit none
  integer, parameter :: rk = 8
  integer, parameter :: ik = 8, nini=21_ik, nfin=20_ik,nfoutx=19_ik,nfoute=18_ik
  integer, parameter :: nfoutt=17_ik, nfoutl=100_ik
  integer(ik), parameter :: ndof=100, ndof2=ndof*2._rk
  integer(ik) :: niter, ii, jj, kk, ll, Ecount
  integer(4) :: nseed
  integer(4), allocatable :: seed(:)
  real(rk), dimension(ndof2) :: x, w, pr !x=(q1,q2,...,qN,p1,p2,...,pN)
  real(rk), dimension(ndof) :: atgc, D, a, dxdt, dpdt, p0, xrange, nscale !at=0, gc=1
  real(rk), dimension(ndof) :: ddxdt, ddpdt
  real(rk), dimension(ndof,ndof) :: Dgc, Dat, agc, aat
  real(rk), dimension(9) :: cc,dc
  real(rk) :: mu, sigma, pi=3.1415926535
  real(rk) :: h, tf, E0, ma, kb, lambda, En, time1, time2, chi, alpha1=0._rk
  character(len=8) :: x1, fmt
  fmt = '(I2.2)'
  ! b same
  ! K = 0.06_rk
  ! rho = 1._rk
  ma = 0.031_rk

  ! D = 0.03_rk
  ! a = 4.5_rk
  ! ma = 300._rk
  kb = 0.00008617_rk
  call cpu_time(time1)
  h = 0.015_rk
  tf = 10000.0_rk
  niter = tf/h
  E0 = 2.0_rk
  ! D and a coefficients
  atgc = 1._rk
  atgc(ndof-4:ndof) = 0._rk
  Dgc = 0._rk
  Dat = 0._rk
  agc = 0._rk
  aat = 0._rk
  do ii = 1,ndof
    Dgc(ii,ii) = 0.075_rk
    Dat(ii,ii) = 0.05_rk
    agc(ii,ii) = 6.9_rk
    aat(ii,ii) = 4.2_rk
  end do
  D = MATMUL(Dgc,atgc) + MATMUL(Dat, 1._rk - atgc)
  a = MATMUL(agc,atgc) + MATMUL(aat, 1._rk - atgc)

  call init_random_seed()
  ! do Ecount=1,8
    ! E0 = Ecount * 1._rk
    ! write (x1,fmt) Ecount
    mu = 50._rk
    sigma = 16._rk
    do ii = 1,ndof
      xrange(ii) = ii*1._rk
    end do
    nscale = 1._rk/(sqrt(2._rk*pi*sigma**2))*exp((-(xrange-mu)**2/sigma**2))
    ! write initial conditions (for posterity)
    open(unit=nini, file='dna_atgc_aba1064.ini')
    write(nini,*) E0
    do jj=1,ndof
      write(nini,*) 0._rk
    end do

    ! for random ICs
    call random_number(pr)
    p0 = sqrt (-2._rk * log(pr(1:ndof))) * cos(2._rk * pi * pr(ndof+1:ndof2))
    ! for SS excitations
    ! p0 = 0._rk
    ! p0(50) = 1._rk
    do jj=1,ndof
      write(nini,*) p0(jj)
    end do

    ! for specific ICs

    ! do jj=1,ndof/2 - 1
    !   write(nini,*) 0
    ! end do
    ! write(nini,*) 1._rk
    ! do jj=ndof/2 +1,ndof
    !   write(nini,*) 0
    ! end do

    ! write(nini,*) 1._rk
    ! write(nini,*) 0._rk
    close(nini)
    open(unit=nfin, file='dna_atgc_aba1064.ini')
    read(nfin,*) E0
    do ii=1,ndof2
      read(nfin,*) x(ii)
    end do
    close(nfin)



    lambda = E0 * 2._rk*ma /(norm2(x(ndof+1:ndof2))**2)
    x(ndof+1:ndof2) = x(ndof+1:ndof2)*sqrt(lambda)

    ! initialise deviation vectors
    call init_random_seed()
    ! call random_number(w)
    ! w(ndof+1:ndof2) = 0._rk
    w = 1._rk
    w = w/(norm2(w))

    ! order 4 coefficients
    ! cc(1) = 1._rk/(2._rk*(2._rk-2._rk**(1._rk/3._rk)))
    ! cc(2) = (1._rk-2._rk**(1._rk/3._rk)) / (2._rk*(2._rk-2._rk**(1._rk/3._rk)))
    ! cc(3) = cc(2)
    ! cc(4) = cc(1)
    !
    ! dc(1) = 1._rk / (2._rk - 2._rk**(1._rk/3._rk))
    ! dc(2) = - (2._rk**(1._rk/3._rk))/(2._rk - 2._rk**(1._rk/3._rk))
    ! dc(3) = dc(1)
    ! dc(4) = 0._rk

    ! aba1064 coefficients
    cc(1) = 0.03809449742241219545697532230863756534060_rk
    cc(2) = 0.1452987161169137492940200726606637497442_rk
    cc(3) = 0.2076276957255412507162056113249882065158_rk
    cc(4) = 0.4359097036515261592231548624010651844006_rk
    cc(5) = -0.6538612258327867093807117373907094120024_rk
    cc(6) = cc(4)
    cc(7) = cc(3)
    cc(8) = cc(2)
    cc(9) = cc(1)

    dc(1) = 0.09585888083707521061077150377145884776921_rk
    dc(2) = 0.2044461531429987806805077839164344779763_rk
    dc(3) = 0.2170703479789911017143385924306336714532_rk
    dc(4) = -0.01737538195906509300561788011852699719871_rk
    dc(5) = dc(4)
    dc(6) = dc(3)
    dc(7) = dc(2)
    dc(8) = dc(1)
    dc(9) = 0._rk
    print*, ' '
    ! open(unit=nfoutx, file='Dna_Order4_xout.csv')
    open(unit=nfoute, file='Dna_ABA1064_Eout.csv')
    open(unit=nfoutt, file='Dna_ABA1064_Tout.csv')
    open(unit=nfoutl, file='Dna_ABA1064_Lout.csv')

    ! open(unit=nfoutt, file='Dna_Order4_Tout'//trim(x1)//'.csv')
    ! open(unit=nfoutx, file='Dna_Order4_xout.csv')

    write(nfoute,*) 0._rk
    write(nfoutt,*) temp(x)
    write(nfoutl,*) 0._rk
    ! write(nfoutx,*) x(1:ndof)
    kk = 0_ik
    do while (kk .lt. tf/h)
      do ll=1,size(cc)
        dxdt = ham_eqs_x(x)
        x(1:ndof) = x(1:ndof) + h * cc(ll) * dxdt(1:ndof)
        ddxdt = var_eqs_x(w)
        w(1:ndof) = w(1:ndof) + h * cc(ll) * ddxdt(1:ndof)
        dpdt = ham_eqs_p(x,D,a)
        x(ndof+1:ndof2) = x(ndof+1:ndof2) - h * dc(ll) * dpdt(1:ndof)
        ddpdt = var_eqs_p(x,w,D,a)
        w(ndof+1:ndof2) = w(ndof+1:ndof2) - h * dc(ll) * ddpdt(1:ndof)
      end do
      ! print*,' x=', w(3), ' p=', w(103), 'lnorm=', log10(norm2(w))
      ! print*, 'dx=', ddxdt(3), 'dp=', ddpdt(3)
      En = energy(x,D,a)
      alpha1 = alpha1 + log(norm2(w))
      chi = (1._rk/(kk*h))*alpha1

      w = w/norm2(w)
      if (mod(kk,10)==0_ik) then

        write(nfoute,*) abs((En-E0)/E0), ', ', maxval(x(1:ndof))
        write(nfoutt,*) temp(x)
        write(nfoutl,*) chi
        ! write(nfoutx,*) x(1:ndof)
      end if
      kk = kk + 1_ik
    end do
    close(nfoute)
    close(nfoutt)
    close(nfoutl)
    ! close(nfoutx)
    call cpu_time(time2)
    print*, (time2-time1)
  ! end do
end program dnao4
