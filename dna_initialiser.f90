program dnao4
  use ham_eqs_mod, only: ham_eqs_x, ham_eqs_p
  use energy_eqns, only: energy, temp
  use var_eqs_mod, only: var_eqs_x, var_eqs_p
  use rand_mod, only: init_random_seed
  use integrators, only: aba864
  implicit none
  integer, parameter :: rk = 8
  integer, parameter :: ik = 8, nini=21_ik, nfin=20_ik
  integer(ik), parameter :: ndof=100, ndof2=ndof*2
  integer(ik) :: ii, jj!,niter
  real(rk), dimension(ndof2) :: x, w, pr !x=(q1,q2,...,qN,p1,p2,...,pN)
  real(rk), dimension(ndof) :: atgc, D, a, p0, xrange !at=0, gc=1
  real(rk), dimension(ndof,ndof) :: Dgc, Dat, agc, aat
  real(rk) :: pi=3.1415926535
  real(rk) :: h, tf, E0, ma, kb, lambda
  ! character(len=8) :: x1, fmt
  ! fmt = '(I2.2)'
  ! b same
  ! K = 0.06_rk
  ! rho = 1._rk
  ma = 0.031_rk

  ! D = 0.03_rk
  ! a = 4.5_rk
  ! ma = 300._rk
  kb = 0.00008617_rk
  h = 0.015_rk
  tf = 1000.0_rk
  ! niter = tf/h
  E0 = 4.0_rk
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
    do ii = 1,ndof
      xrange(ii) = ii*1._rk
    end do
    ! write initial conditions (for posterity)
    open(unit=nini, file='dna_atgc_aba864.ini')
    write(nini,*) E0
    do jj=1,ndof
      write(nini,*) 0._rk
    end do

    ! for random ICs
    call random_number(pr)
    ! Box-Muller method
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
    open(unit=nfin, file='dna_atgc_aba864.ini')
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
    call aba864(x,w,D,a,E0,h,tf)
  end program
