module random
use, non_intrinsic :: kinds, only: i32, i64, sp, dp
use, non_intrinsic :: constants, only: twopi_dp
use, non_intrinsic :: system, only: debug_error_condition
implicit none
private
public :: random_uniform_i32, random_normal_sp, random_normal_multi

    interface random_normal_multi
        module procedure :: random_normal_multi_2d_sp
    end interface random_normal_multi

contains

    impure subroutine random_uniform_i32(x, n, lo, hi)
        integer(kind=i32), intent(in) :: n, lo, hi
        integer(kind=i32), intent(out) :: x(n)
        real(kind=dp) :: work(n), hi_lo
        integer(kind=i64) :: lo_i64
        call debug_error_condition(int(n, kind=i64) > huge(1_i32), &
                                   'RANDOM::RANDOM_UNIFORM_I32 supplied n too large for i32 storage')
        call random_number(work)
        lo_i64 = int(lo, kind=i64)
        hi_lo = real(int(hi, kind=i64) - lo_i64 + 1_i64, kind=dp)
        x = int(int(work*hi_lo, kind=i64) + lo_i64, kind=i32)
    end subroutine random_uniform_i32

    impure subroutine random_normal_sp(x, n, mu, sig)
        integer(kind=i32), intent(in) :: n
        real(kind=sp), intent(out) :: x(n)
        real(kind=sp), intent(in) :: mu, sig
        real(kind=dp) :: x_dp(n), u(n), r, theta, mu_dp, sig_dp
        integer(kind=i32) :: n_2, i
        call debug_error_condition(int(n, kind=i64) > huge(1_i32), &
                                   'RANDOM::RANDOM_NORMAL_SP supplied n too large for i32 storage')
        call random_number(u)
        n_2 = n/2_i32
        mu_dp = real(mu, kind=dp)
        sig_dp = real(sig, kind=dp)
        do concurrent (i=1_i32:n_2)
            r = sqrt(-2.0_dp*log(1.0_dp - u(i)))
            theta = twopi_dp*u(i+n_2)
            x_dp(i) = mu_dp + sig_dp*r*cos(theta)
            x_dp(i+n_2) = mu_dp + sig_dp*r*sin(theta)
        end do
        if (mod(n, 2_i32) /= 0_i32) then
            call random_number(u(1:2))
            x_dp(n) = mu_dp + sig_dp*sqrt(-2.0_dp*log(1.0_dp-u(1)))*cos(twopi_dp*u(2))
        end if
        x = real(x_dp, kind=sp)
    end subroutine random_normal_sp

    impure subroutine random_normal_multi_2d_sp(x, mu, cholesky_factor)
        real(kind=sp), intent(out) :: x(:,:)
        real(kind=sp), intent(in) :: mu(:), cholesky_factor(:,:)
        integer(kind=i32) :: i
        call debug_error_condition((size(x, dim=1, kind=i64) > huge(1_i32)) .or. (size(x, dim=2, kind=i64) > huge(1_i32)), &
                                   'RANDOM::RANDOM_NORMAL_MULTI_2D_SP supplied x too large for i32 storage')
        call debug_error_condition(size(x, dim=1) /= size(mu), &
                                   'RANDOM::RANDOM_NORMAL_MULTI_2D_SP mu does not match x dimensions')
        call debug_error_condition((size(x, dim=1) /= size(cholesky_factor, dim=1)) .or. &
                                   (size(x, dim=1) /= size(cholesky_factor, dim=2)), &
                                   'RANDOM::RANDOM_NORMAL_MULTI_2D_SP Cholesky factor does not match x dimensions')
        call debug_error_condition(size(cholesky_factor, dim=1) /= size(cholesky_factor, dim=2), &
                                   'RANDOM::RANDOM_NORMAL_MULTI_2D_SP Cholesky factor must be square matrix')
        call debug_error_condition(any([(cholesky_factor(i,i), i=1,size(cholesky_factor,dim=1))] < 0.0_sp), &
                                   'RANDOM::RANDOM_NORMAL_MULTI_2D_SP Cholesky factor malformed with negative diagonal')
    end subroutine random_normal_multi_2d_sp

end module random
