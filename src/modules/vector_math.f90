module vector_math
use, non_intrinsic :: kinds, only: i32, i64, sp, dp
use, non_intrinsic :: constants, only: vec_len
use, non_intrinsic :: system, only: debug_error_condition, nearly
implicit none
private

    interface vdot
        module procedure :: vdot_i32
        module procedure :: vdot_i64
        module procedure :: vdot_sp
        module procedure :: vdot_dp
    end interface vdot

    interface vmag
        module procedure :: vmag_sp
        module procedure :: vmag_dp
    end interface vmag

    interface vunit
        module procedure :: vunit_sp
        module procedure :: vunit_inplace_sp
        module procedure :: vunit_dp
        module procedure :: vunit_inplace_dp
    end interface vunit

    public :: vdot, vmag, vunit

contains

    pure function vdot_i32(v1, v2) result(val)
        integer(kind=i32), intent(in) :: v1(:), v2(size(v1, kind=i64))
        integer(kind=i32) :: val
        integer(kind=i32) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 0_i32
        n = size(v1, kind=i64)
        if (n > 0_i64) then
            accumulators = 0_i32
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators + v1(i:i+vec_len-1_i64)*v2(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + v1(n-i+1_i64:n)*v2(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function vdot_i32

    pure function vdot_i64(v1, v2) result(val)
        integer(kind=i64), intent(in) :: v1(:), v2(size(v1, kind=i64))
        integer(kind=i64) :: val
        integer(kind=i64) :: n, i, accumulators(vec_len)
        val = 0_i64
        n = size(v1, kind=i64)
        if (n > 0_i64) then
            accumulators = 0_i64
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators + v1(i:i+vec_len-1_i64)*v2(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + v1(n-i+1_i64:n)*v2(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function vdot_i64

    pure function vdot_sp(v1, v2) result(val)
        real(kind=sp), intent(in) :: v1(:), v2(size(v1, kind=i64))
        real(kind=sp) :: val
        real(kind=sp) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 0.0_sp
        n = size(v1, kind=i64)
        if (n > 0_i64) then
            accumulators = 0.0_sp
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators + v1(i:i+vec_len-1_i64)*v2(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + v1(n-i+1_i64:n)*v2(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function vdot_sp

    pure function vdot_dp(v1, v2) result(val)
        real(kind=dp), intent(in) :: v1(:), v2(size(v1, kind=i64))
        real(kind=dp) :: val
        real(kind=dp) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 0.0_dp
        n = size(v1, kind=i64)
        if (n > 0_i64) then
            accumulators = 0.0_dp
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators + v1(i:i+vec_len-1_i64)*v2(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + v1(n-i+1_i64:n)*v2(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function vdot_dp

    pure function vmag_sp(v) result(val)
        real(kind=sp), intent(in) :: v(:)
        real(kind=sp) :: val
        val = sqrt(vdot(v,v))
    end function vmag_sp

    pure function vmag_dp(v) result(val)
        real(kind=dp), intent(in) :: v(:)
        real(kind=dp) :: val
        val = sqrt(vdot(v,v))
    end function vmag_dp

    pure subroutine vunit_sp(v, vhat)
        real(kind=sp), intent(in) :: v(:)
        real(kind=sp), intent(out) :: vhat(:)
        real(kind=sp) :: mag
        mag = vmag(v)
        call debug_error_condition(nearly(mag, 0.0_sp), &
                                   'module VECTOR_MATH :: vunit subroutine invalid for vectors with magnitude near 0.0')
        vhat = v/mag
    end subroutine vunit_sp

    pure subroutine vunit_inplace_sp(v)
        real(kind=sp), intent(inout) :: v(:)
        real(kind=sp) :: mag
        mag = vmag(v)
        call debug_error_condition(nearly(mag, 0.0_sp), &
                                   'module VECTOR_MATH :: vunit subroutine invalid for vectors with magnitude near 0.0')
        v = v/mag
    end subroutine vunit_inplace_sp

    pure subroutine vunit_dp(v, vhat)
        real(kind=dp), intent(in) :: v(:)
        real(kind=dp), intent(out) :: vhat(:)
        real(kind=dp) :: mag
        mag = vmag(v)
        call debug_error_condition(nearly(mag, 0.0_dp), &
                                   'module VECTOR_MATH :: vunit subroutine invalid for vectors with magnitude near 0.0')
        vhat = v/mag
    end subroutine vunit_dp

    pure subroutine vunit_inplace_dp(v)
        real(kind=dp), intent(inout) :: v(:)
        real(kind=dp) :: mag
        mag = vmag(v)
        call debug_error_condition(nearly(mag, 0.0_dp), &
                                   'module VECTOR_MATH :: vunit subroutine invalid for vectors with magnitude near 0.0')
        v = v/mag
    end subroutine vunit_inplace_dp

end module vector_math
