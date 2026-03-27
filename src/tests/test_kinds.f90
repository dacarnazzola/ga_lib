program test_kinds
use, non_intrinsic :: kinds, only: stdin, stdout, stderr, i32, i64, sp, dp, c_bool, c_char
implicit none

    character(len=*, kind=c_char), parameter :: fmt1 = '(a,i0)'
    character(len=*, kind=c_char), parameter :: fmt2 = '(a,i0,a,i0,a)'

    integer(kind=i32) :: test_i32
    integer(kind=i64) :: test_i64
    real(kind=sp) :: test_sp
    real(kind=dp) :: test_dp
    logical(kind=c_bool) :: test_c_bool
    character(kind=c_char) :: test_c_char

    write(unit=stdout, fmt=fmt1) 'TEST_KINDS :: stdin = ',stdin
    write(unit=stdout, fmt=fmt1) 'TEST_KINDS :: stdout = ',stdout
    write(unit=stdout, fmt=fmt1) 'TEST_KINDS :: stderr = ',stderr

    write(unit=stdout, fmt=fmt2) 'TEST_KINDS :: integer(kind=i32) = integer(kind=',i32,'), ',storage_size(test_i32),' bits'
    write(unit=stdout, fmt=fmt2) 'TEST_KINDS :: integer(kind=i64) = integer(kind=',i64,'), ',storage_size(test_i64),' bits'
    write(unit=stdout, fmt=fmt2) 'TEST_KINDS :: real(kind=i32) = real(kind=',sp,'), ',storage_size(test_sp),' bits'
    write(unit=stdout, fmt=fmt2) 'TEST_KINDS :: real(kind=i32) = real(kind=',dp,'), ',storage_size(test_dp),' bits'
    write(unit=stdout, fmt=fmt2) 'TEST_KINDS :: logical(kind=c_bool)  = logical(kind=',c_bool,'), ', &
                                 storage_size(test_c_bool),' bits'
    write(unit=stdout, fmt=fmt2) 'TEST_KINDS :: character(kind=c_char) = character(kind=',c_char,'), ', &
                                 storage_size(test_c_char),' bits'

end program test_kinds
