program test_random
use, non_intrinsic :: kinds, only: stdout, i32, i64, dp
use, non_intrinsic :: random, only: random_uniform_i32
implicit none

    integer(kind=i32), parameter :: n = 100000, xlo = 1, xhi = 6

    integer(kind=i32) :: x(n), i

    call random_uniform_i32(x, size(x, kind=i32), xlo, xhi)

    !$omp parallel do
    do i=xlo,xhi
        write(unit=stdout, fmt='(a,i0,a,f0.1,a)') 'x==',i,' (',count(x==i, kind=i64)*100.0_dp/size(x, kind=i64),'%)'
    end do
    !$omp end parallel do

end program test_random
