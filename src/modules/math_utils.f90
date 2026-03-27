module math_utils
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, c_bool
use, non_intrinsic :: system, only: debug_error_condition
implicit none
private

    interface clamp
        module procedure :: clamp_i32
        module procedure :: clamp_i64
        module procedure :: clamp_sp
        module procedure :: clamp_dp
    end interface clamp

    interface linspace
        module procedure :: linspace_sp
        module procedure :: linspace_dp
    end interface linspace

    public :: clamp, linspace

contains

    pure elemental subroutine clamp_i32(x, xlo, xhi)
        integer(kind=i32), intent(inout) :: x
        integer(kind=i32), intent(in) :: xlo, xhi
        x = min(max(x, xlo), xhi)
    end subroutine clamp_i32

    pure elemental subroutine clamp_i64(x, xlo, xhi)
        integer(kind=i64), intent(inout) :: x
        integer(kind=i64), intent(in) :: xlo, xhi
        x = min(max(x, xlo), xhi)
    end subroutine clamp_i64

    pure elemental subroutine clamp_sp(x, xlo, xhi)
        real(kind=sp), intent(inout) :: x
        real(kind=sp), intent(in) :: xlo, xhi
        x = min(max(x, xlo), xhi)
    end subroutine clamp_sp

    pure elemental subroutine clamp_dp(x, xlo, xhi)
        real(kind=dp), intent(inout) :: x
        real(kind=dp), intent(in) :: xlo, xhi
        x = min(max(x, xlo), xhi)
    end subroutine clamp_dp

    pure subroutine linspace_sp(x, xmin, xmax)
        real(kind=sp), intent(out) :: x(:)
        real(kind=sp), intent(in) :: xmin, xmax
        integer(kind=i64) :: n, i
        real(kind=sp) :: dx
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module MATH_UTILS :: subroutine linspace requires n >= 1')
        dx = (xmax - xmin)/real(n - 1_i64, kind=sp)
        x(1_i64) = xmin
        do i=2_i64,n-1_i64
            x(i) = xmin + dx*real(i - 1_i64, kind=sp)
        end do
        x(n) = xmax
    end subroutine linspace_sp

    pure subroutine linspace_dp(x, xmin, xmax)
        real(kind=dp), intent(out) :: x(:)
        real(kind=dp), intent(in) :: xmin, xmax
        integer(kind=i64) :: n, i
        real(kind=dp) :: dx
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module MATH_UTILS :: subroutine linspace requires n >= 1')
        dx = (xmax - xmin)/real(n - 1_i64, kind=dp)
        x(1_i64) = xmin
        do i=2_i64,n-1_i64
            x(i) = xmin + dx*real(i - 1_i64, kind=dp)
        end do
        x(n) = xmax
    end subroutine linspace_dp

end module math_utils
