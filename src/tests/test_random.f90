program test_random
use, non_intrinsic :: kinds, only: stdout, i32, i64, sp, dp
use, non_intrinsic :: random, only: random_uniform, random_normal
use, non_intrinsic :: statistics, only: avg, std
implicit none

    integer(kind=i32), allocatable :: xi32(:)
    integer(kind=i64), allocatable :: xi64(:)
    real(kind=sp), allocatable :: xsp(:)
    real(kind=dp), allocatable :: xdp(:)
    integer(kind=i32) :: i

    do i=0_i32,3_i32

        if (allocated(xi32)) deallocate(xi32); allocate(xi32(i))
        if (allocated(xi64)) deallocate(xi64); allocate(xi64(i))
        if (allocated(xsp)) deallocate(xsp); allocate(xsp(i))
        if (allocated(xdp)) deallocate(xdp); allocate(xdp(i))

        write(unit=stdout, fmt='(a,i0,a)') 'TEST_RANDOM :: arrays with ',i,' elements'

        call random_uniform(xi32, 1_i32, 6_i32)
        write(unit=stdout, fmt='(a,i0," ",i0)') 'TEST_RANDOM :: random_uniform xi32 min/max ',minval(xi32),maxval(xi32)

        call random_uniform(xi64, 1_i64, 6_i64)
        write(unit=stdout, fmt='(a,i0," ",i0)') 'TEST_RANDOM :: random_uniform xi64 min/max ',minval(xi64),maxval(xi64)

        call random_uniform(xsp, 1.0_sp, 6.0_sp)
        write(unit=stdout, fmt='(a,f0.1," ",f0.1)') 'TEST_RANDOM :: random_uniform xsp min/max ',minval(xsp),maxval(xsp)

        call random_uniform(xdp, 1.0_dp, 6.0_dp)
        write(unit=stdout, fmt='(a,f0.1," ",f0.1)') 'TEST_RANDOM :: random_uniform xdp min/max ',minval(xdp),maxval(xdp)

        call random_normal(xsp, 100.0_sp, 10.0_sp)
        write(unit=stdout, fmt='(a,f0.1," ",f0.1)') 'TEST_RANDOM :: random_normal xsp min/max ',minval(xsp),maxval(xsp)

        call random_normal(xdp, 100.0_dp, 10.0_dp)
        write(unit=stdout, fmt='(a,f0.1," ",f0.1,/)') 'TEST_RANDOM :: random_normal xdp min/max ',minval(xdp),maxval(xdp)

    end do

    i = 2_i32**20_i32
    if (allocated(xsp)) deallocate(xsp); allocate(xsp(i))
    if (allocated(xdp)) deallocate(xdp); allocate(xdp(i))

    write(unit=stdout, fmt='(a,i0,a)') 'TEST_RANDOM :: arrays with ',i,' elements'

    call random_normal(xsp, 19.93_sp, 8.31_sp)
    write(unit=stdout, fmt='(a,f0.2," ",f0.2)') 'TEST_RANDOM :: random_normal xsp avg/std ',avg(xsp),std(xsp)

    call random_normal(xdp, 19.94_dp, 5.30_dp)
    write(unit=stdout, fmt='(a,f0.2," ",f0.2)') 'TEST_RANDOM :: random_normal xdp avg/std ',avg(xdp),std(xdp)

end program test_random
