module integrators
  implicit none
  public
    ! aba864 coefficients
    cc(1) = 0.0711334264982231177779387300061549964174_rk
    cc(2) = 0.241153427956640098736487795326289649618_rk
    cc(3) = 0.521411761772814789212136078067994229991_rk
    cc(4) = -0.333698616227678005726562603400438876027_rk
    cc(5) = cc(4)
    cc(6) = cc(3)
    cc(7) = cc(2)
    cc(8) = cc(1)

    dc(1) = 0.183083687472197221961703757166430291072_rk
    dc(2) = 0.310782859898574869507522291054262796375_rk
    dc(3) = -0.0265646185119588006972121379164987592663_rk
    dc(4) = 0.0653961422823734184559721793911134363710_rk
    dc(5) = dc(3)
    dc(6) = dc(2)
    dc(7) = dc(1)
    dc(8) = 0._rk
    print*, ' '
    ! open(unit=nfoutx, file='Dna_Order4_xout.csv')
    open(unit=nfoute, file='Dna_ABA864_Eout.csv')
    open(unit=nfoutt, file='Dna_ABA864_Tout.csv')
    open(unit=nfoutl, file='Dna_ABA864_Lout.csv')

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
