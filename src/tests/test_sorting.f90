program test_sorting
use, non_intrinsic :: kinds, only: stdout, i32, i64, sp, dp
use, non_intrinsic :: random, only: random_uniform
use, non_intrinsic :: sorting, only: is_sorted, sort
use, non_intrinsic :: system, only: nearly
implicit none

    integer(kind=i64), parameter :: n = 128_i64

    integer(kind=i32) :: xi32(n), xi322(n)
    integer(kind=i64) :: xi64(n), xindi64(n), i, xi642(n)
    real(kind=sp) :: xsp(n), xsp2(n)
    real(kind=dp) :: xdp(n), xdp2(n)

    call random_uniform(xi32, -100_i32, 100_i32)
    write(unit=stdout, fmt='(a,l1)') 'TEST_SORTING :: random i32 array is_sorted: ',is_sorted(xi32)
    write(unit=stdout, fmt='(a,*(i0," "))') 'TEST_SORTING :: i32 array: ',xi32
    xi322 = xi32
    call sort(xi32)
    write(unit=stdout, fmt='(a,l1)') 'TEST_SORTING :: i32 array after SORT is_sorted: ',is_sorted(xi32)
    write(unit=stdout, fmt='(a,*(i0," "))') 'TEST_SORTING :: i32 array: ',xi32
    xi32 = xi322
    xindi64 = [(i, i=1_i64,n)]
    call sort(xi32, xindi64)
    write(unit=stdout, fmt='(2(a,l1),/)') 'TEST_SORTING :: i32 array after SORT with indices is_sorted: ',is_sorted(xi32), &
                                       ' and sorted indices: ',all(xi32 == xi322(xindi64))

    call random_uniform(xi64, -100_i64, 100_i64)
    write(unit=stdout, fmt='(a,l1)') 'TEST_SORTING :: random i64 array is_sorted: ',is_sorted(xi64)
    write(unit=stdout, fmt='(a,*(i0," "))') 'TEST_SORTING :: i64 array: ',xi64
    xi642 = xi64
    call sort(xi64)
    write(unit=stdout, fmt='(a,l1)') 'TEST_SORTING :: i64 array after SORT is_sorted: ',is_sorted(xi64)
    write(unit=stdout, fmt='(a,*(i0," "))') 'TEST_SORTING :: i64 array: ',xi64
    xi64 = xi642
    xindi64 = [(i, i=1_i64,n)]
    call sort(xi64, xindi64)
    write(unit=stdout, fmt='(2(a,l1),/)') 'TEST_SORTING :: i64 array after SORT with indices is_sorted: ',is_sorted(xi64), &
                                       ' and sorted indices: ',all(xi64 == xi642(xindi64))

    call random_uniform(xsp, -100.0_sp, 100.0_sp)
    write(unit=stdout, fmt='(a,l1)') 'TEST_SORTING :: random sp array is_sorted: ',is_sorted(xsp)
    write(unit=stdout, fmt='(a,*(f0.1," "))') 'TEST_SORTING :: sp array: ',xsp
    xsp2 = xsp
    call sort(xsp)
    write(unit=stdout, fmt='(a,l1)') 'TEST_SORTING :: sp array after SORT is_sorted: ',is_sorted(xsp)
    write(unit=stdout, fmt='(a,*(f0.1," "))') 'TEST_SORTING :: sp array: ',xsp
    xsp = xsp2
    xindi64 = [(i, i=1_i64,n)]
    call sort(xsp, xindi64)
    write(unit=stdout, fmt='(2(a,l1),/)') 'TEST_SORTING :: sp array after SORT with indices is_sorted: ',is_sorted(xsp), &
                                       ' and sorted indices: ',all(nearly(xsp,xsp2(xindi64)))

    call random_uniform(xdp, -100.0_dp, 100.0_dp)
    write(unit=stdout, fmt='(a,l1)') 'TEST_SORTING :: random dp array is_sorted: ',is_sorted(xdp)
    write(unit=stdout, fmt='(a,*(f0.1," "))') 'TEST_SORTING :: dp array: ',xdp
    xdp2 = xdp
    call sort(xdp)
    write(unit=stdout, fmt='(a,l1)') 'TEST_SORTING :: dp array after SORT is_sorted: ',is_sorted(xdp)
    write(unit=stdout, fmt='(a,*(f0.1," "))') 'TEST_SORTING :: dp array: ',xdp
    xdp = xdp2
    xindi64 = [(i, i=1_i64,n)]
    call sort(xdp, xindi64)
    write(unit=stdout, fmt='(2(a,l1),/)') 'TEST_SORTING :: dp array after SORT with indices is_sorted: ',is_sorted(xdp), &
                                       ' and sorted indices: ',all(nearly(xdp,xdp2(xindi64)))

end program test_sorting
