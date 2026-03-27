program test_debug
use, non_intrinsic :: kinds, only: stdout, dp
use, non_intrinsic :: constants, only: debug
use, non_intrinsic :: timing, only: timer, tic, toc, get_elapsed
implicit none

    type(timer) :: bench
    real(kind=dp) :: elapsed

    write(unit=stdout, fmt='(a,l1)') 'TEST_DEBUG :: status of logical parameter DEBUG = ',debug

    write(unit=stdout, fmt='(a)') 'TEST_DEBUG :: about to call toc(timer) before tic(timer) --> expect error stop if DEBUG=T'
    call toc(bench)
    
    write(unit=stdout, fmt='(a)') 'TEST_DEBUG :: about to call get_elapsed(timer) before tic(timer) or toc(timer) '// &
                                  '--> expect error stop if DEBUG=T'
    elapsed = get_elapsed(bench)
    write(unit=stdout, fmt='(a,g0)') 'TEST_DEBUG :: value of elapsed should be NaN or Infinity, elapsed=',elapsed

end program test_debug
