module kinds
use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit
use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_float, c_double, c_bool, c_char
implicit none
private

    integer, parameter :: stdin  = input_unit
    integer, parameter :: stdout = output_unit
    integer, parameter :: stderr = error_unit

    integer, parameter :: i32 = c_int32_t
    integer, parameter :: i64 = c_int64_t
    integer, parameter :: sp  = c_float
    integer, parameter :: dp  = c_double

    public :: stdin, stdout, stderr, i32, i64, sp, dp, c_bool, c_char

end module kinds
