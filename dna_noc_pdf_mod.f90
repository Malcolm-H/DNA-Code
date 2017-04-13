module noc_pdf_mod
  use rand_mod, only: init_random_seed
  use noc_distribution_mod, only: get_atgc
  implicit none
CONTAINS
  subroutine write_pdf_data(ndof, pat, niter, fmt2)
    implicit none
    integer, parameter :: rk=8, ik=4
    integer(ik) :: ii, jj, kk, ndof, niter, acount, flag, noc
    real(rk) :: pat, t1, t2
    integer(ik), dimension(ndof) :: atgc
    character(len=2*ndof) :: atgcchar
    character(len=ndof) :: at1, at2
    character(len=8) :: x1, x2, fmt, fmt2

    fmt = '(I3)'
    write(x1,fmt2) niter
    write(x2,fmt) ndof
    call init_random_seed()
    !!! choose atgc randomly based on ndof, ngc
    open(unit=122, file='../DNA-Data/ATGC_'//trim(x2)//'_n'//trim(x1)//'.txt')
    ! open(unit=121, file='../DNA-Data/time_noc.txt')
    ! do kk=50,3000,50
      acount = 0
      open(unit=123, status='scratch', access='direct', recl=200, action='readwrite')
      ! call cpu_time(t1)
      do jj=1,niter
        flag = 0
        atgc = get_atgc(ndof, (ndof- int(pat * ndof)))
        noc = 0_ik
        do ii=1,ndof-1
          if(atgc(ii) .ne. atgc(ii+1)) noc = noc + 1_ik
        end do
        if (atgc(1) .ne. atgc(ndof)) noc = noc + 1_ik
        ! this and the two loops create the concatenated string to be checked against
        !-----------------------------------------------------
        atgcchar = char(int(atgc(1)+48))
        do ii=2,10
          atgcchar = trim(atgcchar) // trim(char(int(atgc(ii) + 48)))
        end do
        do ii=1,10
          atgcchar = trim(atgcchar) // trim(char(int(atgc(ii) + 48)))
        end do
        !-----------------------------------------------------
        ! loop and read recorded atgc values to compare against
        at2 = char(int(atgc(1) + 48))
        do ii=2,ndof
          at2 = trim(at2) // trim(char(int(atgc(ii) + 48)))
        end do
        if (jj .gt. 1) then
          do ii = 1,acount
            read(unit=123,rec=ii) at1
            atgcchar = trim(at1) // trim(at1)

            if (index(atgcchar, at2) .ne. 0) then
              flag = 1
              print*, 1
              exit
            end if
          end do
        end if
        !-----------------------------------------------------
        if (flag .eq. 0) then
          write(unit=123,rec=jj) at2
          write(122,*) at2
          acount = acount + 1
        end if
        atgc = 0._rk
      end do
      ! call cpu_time(t2)
      ! write(121,*) kk, ', ', (t2-t1)
      close(123)

    ! end do
    ! close(121)
    close(122)
  end subroutine
end module
