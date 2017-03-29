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
