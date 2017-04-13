module noc_distribution_mod
  use rand_mod, only: init_random_seed
  implicit none
  public :: noc_distribution, get_atgc
CONTAINS
  function noc_distribution(noc, pat, ntot) result(atgc)
    use rand_mod, only: init_random_seed
    implicit none
    integer, parameter :: ik=4, rk=8
    integer(ik), intent(in) :: noc, ntot
    real(rk), intent(in) :: pat
    integer(ik), dimension(ntot) :: atgc
    integer(ik) :: nat, ngc
    integer(ik), dimension(noc/2) :: part_at_size, part_gc_size
    integer(ik) :: ii, jj, kk=1
    real(rk) :: at_size, gc_size
    integer(ik) :: at_int_size, gc_int_size
    nat = int(pat * ntot)
    ngc = int((1-pat) * ntot)
    part_at_size = 1_ik
    part_gc_size = 1_ik
    call init_random_seed()
    do ii=1,(nat-noc/2)
      call random_number(at_size)
      at_int_size = int(at_size*(noc/2)) + 1
      part_at_size(at_int_size) = part_at_size(at_int_size) + 1_ik
    end do
    do ii=1,(ngc-noc/2)
      call random_number(gc_size)
      gc_int_size = int(gc_size*(noc/2)) + 1
      part_gc_size(gc_int_size) = part_gc_size(gc_int_size) + 1_ik
    end do
    ! print*, (sum(part_gc_size)+sum(part_at_size))
    do ii=1,noc
      if(mod(ii,2) .eq. 1_ik) then
        do jj=1,part_at_size((ii+1)/2)
          atgc(kk) = 1_ik
          kk = kk + 1_ik
        end do
      else
        do jj=1,part_gc_size(ii/2)
          atgc(kk) = 0_ik
          kk = kk + 1_ik
        end do
      end if
    end do
  end function noc_distribution
  function get_atgc(ndof, ngc) result(atgc)
    implicit none
    integer, parameter :: rk=8, ik=4
    integer(ik) :: ndof, ngc, rand_ind_i
    integer(ik), dimension(ndof) :: atgc
    real(rk) :: rand_ind_r
    ! call init_random_seed()
    atgc = 0_ik
    do while (sum(atgc) .lt. ngc)
      call random_number(rand_ind_r)
      rand_ind_i = int(ndof * rand_ind_r) + 1_ik
      ! if (atgc(rand_ind_i) .ne. 1)
      atgc(rand_ind_i) = 1_ik
    end do
  end function get_atgc
end module
