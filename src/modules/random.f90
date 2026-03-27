module random
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, c_bool
use, non_intrinsic :: constants, only: twopi_sp, twopi_dp
use, non_intrinsic :: system, only: debug_error_condition
use, non_intrinsic :: matrix_math, only: chol
implicit none
private

    interface random_uniform
        module procedure :: random_uniform_i32
        module procedure :: random_uniform_i64
        module procedure :: random_uniform_sp
        module procedure :: random_uniform_dp
    end interface random_uniform

    interface random_normal
        module procedure :: random_normal_sp
        module procedure :: random_normal_dp
        module procedure :: random_multivariable_normal_sp
        module procedure :: random_multivariable_normal_dp
    end interface random_normal

    interface random_log_uniform
        module procedure :: random_log_uniform_sp
        module procedure :: random_log_uniform_dp
    end interface random_log_uniform

    public :: random_uniform, random_normal, random_log_uniform

contains

    impure subroutine random_uniform_i32(vout, vmin, vmax)
        integer(kind=i32), intent(out) :: vout(:)
        integer(kind=i32), intent(in) :: vmin, vmax
        real(kind=sp) :: work(size(vout, kind=i64))
        call random_number(work)
        vout = int(work*real(vmax - vmin + 1_i64, kind=sp), kind=i32) + vmin
    end subroutine random_uniform_i32

    impure subroutine random_uniform_i64(vout, vmin, vmax)
        integer(kind=i64), intent(out) :: vout(:)
        integer(kind=i64), intent(in) :: vmin, vmax
        real(kind=sp) :: work(size(vout, kind=i64))
        call random_number(work)
        vout = int(work*real(vmax - vmin + 1_i64, kind=sp), kind=i64) + vmin
    end subroutine random_uniform_i64

    impure subroutine random_uniform_sp(vout, vmin, vmax)
        real(kind=sp), intent(out) :: vout(:)
        real(kind=sp), intent(in) :: vmin, vmax
        call random_number(vout)
        vout = vout*(vmax - vmin) + vmin
    end subroutine random_uniform_sp

    impure subroutine random_uniform_dp(vout, vmin, vmax)
        real(kind=dp), intent(out) :: vout(:)
        real(kind=dp), intent(in) :: vmin, vmax
        call random_number(vout)
        vout = vout*(vmax - vmin) + vmin
    end subroutine random_uniform_dp

    impure subroutine random_normal_sp(vout, mu, sig)
        real(kind=sp), intent(out) :: vout(:)
        real(kind=sp), intent(in) :: mu, sig
        real(kind=sp) :: u((size(vout,kind=i64)+1_i64)/2_i64), r((size(vout,kind=i64)+1_i64)/2_i64)
        integer(kind=i64) :: n
        call random_number(r)
        r = 1.0_sp - r
        r = sqrt(-2.0_sp*log(r))
        call random_number(u)
        n = size(vout, kind=i64)
        vout(1_i64:(n/2_i64)) = mu + sig*r(1_i64:(n/2_i64))*cos(twopi_sp*u(1_i64:(n/2_i64)))
        vout((n/2_i64+1):n) = mu + sig*r*sin(twopi_sp*u)
    end subroutine random_normal_sp

    impure subroutine random_normal_dp(vout, mu, sig)
        real(kind=dp), intent(out) :: vout(:)
        real(kind=dp), intent(in) :: mu, sig
        real(kind=dp) :: u((size(vout,kind=i64)+1_i64)/2_i64), r((size(vout,kind=i64)+1_i64)/2_i64)
        integer(kind=i64) :: n
        call random_number(r)
        r = 1.0_dp - r
        r = sqrt(-2.0_dp*log(r))
        call random_number(u)
        n = size(vout, kind=i64)
        vout(1_i64:(n/2_i64)) = mu + sig*r(1_i64:(n/2_i64))*cos(twopi_dp*u(1_i64:(n/2_i64)))
        vout((n/2_i64+1):n) = mu + sig*r*sin(twopi_dp*u)
    end subroutine random_normal_dp

    impure subroutine random_multivariable_normal_sp(vout, mu, cov)
        real(kind=sp), intent(out) :: vout(:,:)
        real(kind=sp), intent(in) :: mu(:), cov(:,:)
        integer(kind=i64) :: ndims, nvars, i
        real(kind=sp) :: cov_L(size(cov, dim=1, kind=i64),size(cov, dim=2, kind=i64)), &
                         z(size(vout, kind=i64))
        ndims = size(mu, kind=i64)
        call debug_error_condition(logical(size(vout, dim=1, kind=i64) /= ndims, kind=c_bool), &
                                   'module RANDOM :: subroutine random_normal has mismatched ndims between vout and mu')
        call debug_error_condition(logical(size(cov, dim=1, kind=i64) /= ndims, kind=c_bool), &
                                   'module RANDOM :: subroutine random_normal has mismatched ndims between cov and mu')
        call debug_error_condition(logical(size(cov, dim=1, kind=i64) /= size(cov, dim=2, kind=i64), kind=c_bool), &
                                   'module RANDOM :: subroutine random_normal requires square covariance input (n x n)')
        nvars = size(vout, dim=2, kind=i64)
        call chol(cov, cov_L)
        call random_normal(z, 0.0_sp, 1.0_sp)
        do i=1_i64,nvars
            vout(:,i) = mu + matmul(cov_L, z((i-1_i64)*ndims+1_i64:i*ndims))
        end do
    end subroutine random_multivariable_normal_sp

    impure subroutine random_multivariable_normal_dp(vout, mu, cov)
        real(kind=dp), intent(out) :: vout(:,:)
        real(kind=dp), intent(in) :: mu(:), cov(:,:)
        integer(kind=i64) :: ndims, nvars, i
        real(kind=dp) :: cov_L(size(cov, dim=1, kind=i64),size(cov, dim=2, kind=i64)), &
                         z(size(vout, kind=i64))
        ndims = size(mu, kind=i64)
        call debug_error_condition(logical(size(vout, dim=1, kind=i64) /= ndims, kind=c_bool), &
                                   'module RANDOM :: subroutine random_normal has mismatched ndims between vout and mu')
        call debug_error_condition(logical(size(cov, dim=1, kind=i64) /= ndims, kind=c_bool), &
                                   'module RANDOM :: subroutine random_normal has mismatched ndims between cov and mu')
        call debug_error_condition(logical(size(cov, dim=1, kind=i64) /= size(cov, dim=2, kind=i64), kind=c_bool), &
                                   'module RANDOM :: subroutine random_normal requires square covariance input (n x n)')
        nvars = size(vout, dim=2, kind=i64)
        call chol(cov, cov_L)
        call random_normal(z, 0.0_dp, 1.0_dp)
        do i=1_i64,nvars
            vout(:,i) = mu + matmul(cov_L, z((i-1_i64)*ndims+1_i64:i*ndims))
        end do
    end subroutine random_multivariable_normal_dp

    impure subroutine random_log_uniform_sp(vout, vmin, vmax)
        real(kind=sp), intent(out) :: vout(:)
        real(kind=sp), intent(in) :: vmin, vmax
        call debug_error_condition(logical((vmin <= 0.0_sp) .or. (vmax <= 0.0_sp), kind=c_bool), &
                                   'module RANDOM :: subroutine random_log_uniform invalid for vmin<=0 or vmax<=0')
        call random_number(vout)
        vout = exp(vout*(log(vmax) - log(vmin)) + log(vmin))
    end subroutine random_log_uniform_sp

    impure subroutine random_log_uniform_dp(vout, vmin, vmax)
        real(kind=dp), intent(out) :: vout(:)
        real(kind=dp), intent(in) :: vmin, vmax
        call debug_error_condition(logical((vmin <= 0.0_dp) .or. (vmax <= 0.0_dp), kind=c_bool), &
                                   'module RANDOM :: subroutine random_log_uniform invalid for vmin<=0 or vmax<=0')
        call random_number(vout)
        vout = exp(vout*(log(vmax) - log(vmin)) + log(vmin))
    end subroutine random_log_uniform_dp

end module random
