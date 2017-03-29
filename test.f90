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


program test_rand
  use rand_mod, only: init_random_seed
  implicit none
  integer(4), allocatable :: seed(:)
  integer(4) :: n
  real(8) :: r(5)
  call init_random_seed()
  ! call random_seed(size=n)
  ! allocate(seed(n))
  ! call random_seed(get=seed)
  ! ! print*, seed
  call random_number(r)
  print*,r
  ! seed = 2
  ! call random_seed(put=seed)
  ! call random_number(r)
  print*,r
end program test_rand
