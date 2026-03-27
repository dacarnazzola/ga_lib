module statistics
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, c_bool
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: constants, only: vec_len, eps_sp, eps_dp, twopi_sp, twopi_dp
use, non_intrinsic :: vector_math, only: vdot
implicit none
private

    interface dsum
        module procedure :: dsum_i32
        module procedure :: dsum_i64
        module procedure :: dsum_sp
        module procedure :: dsum_dp
    end interface dsum

    interface prod
        module procedure :: prod_i32
        module procedure :: prod_i64
        module procedure :: prod_sp
        module procedure :: prod_dp
    end interface prod

    interface avg
        module procedure :: avg_sp
        module procedure :: avg_dp
    end interface avg

    interface geomean
        module procedure :: geomean_sp
        module procedure :: geomean_dp
    end interface geomean

    interface std
        module procedure :: std_sp
        module procedure :: std_dp
        module procedure :: std_weighted_sp
        module procedure :: std_weighted_dp
    end interface std

    interface normal_pdf
        module procedure :: normal_pdf_1d_sp
        module procedure :: normal_pdf_1d_dp
    end interface normal_pdf

    interface mv_normal_pdf
        module procedure :: normal_pdf_diagonal_covariance_sp
        module procedure :: normal_pdf_diagonal_covariance_dp
    end interface mv_normal_pdf

    interface normalize
        module procedure :: normalize_inplace_sp
        module procedure :: normalize_sp
        module procedure :: normalize_inplace_dp
        module procedure :: normalize_dp
    end interface normalize

    interface normalize_min_max
        module procedure :: normalize_min_max_inplace_sp
        module procedure :: normalize_min_max_sp
        module procedure :: normalize_min_max_inplace_dp
        module procedure :: normalize_min_max_dp
    end interface normalize_min_max

    interface cumsum
        module procedure :: cumsum_i32
        module procedure :: cumsum_inplace_i32
        module procedure :: cumsum_i64
        module procedure :: cumsum_inplace_i64
        module procedure :: cumsum_sp
        module procedure :: cumsum_inplace_sp
        module procedure :: cumsum_dp
        module procedure :: cumsum_inplace_dp
    end interface cumsum

    interface cov
        module procedure :: cov_sp
        module procedure :: cov_dp
        module procedure :: cov_weighted_sp
        module procedure :: cov_weighted_dp
    end interface cov

    public :: dsum, avg, geomean, std, normal_pdf, normalize, cumsum, prod, mv_normal_pdf, cov, normalize_min_max

contains

    pure function dsum_i32(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val
        integer(kind=i32) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 0_i32
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 0_i32
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators + x(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + x(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function dsum_i32

    pure function dsum_i64(x) result(val)
        integer(kind=i64), intent(in) :: x(:)
        integer(kind=i64) :: val
        integer(kind=i64) :: n, i, accumulators(vec_len)
        val = 0_i64
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 0_i64
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators + x(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + x(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function dsum_i64

    pure function dsum_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp) :: val
        real(kind=sp) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 0.0_sp
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 0.0_sp
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators + x(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + x(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function dsum_sp

    pure function dsum_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: val
        real(kind=dp) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 0.0_dp
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 0.0_dp
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators + x(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + x(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function dsum_dp

    pure function prod_i32(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val
        integer(kind=i32) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 1_i32
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 1_i32
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators * x(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) * x(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val * accumulators(i) * accumulators(i + 1_i64)
            end do
        end if
    end function prod_i32

    pure function prod_i64(x) result(val)
        integer(kind=i64), intent(in) :: x(:)
        integer(kind=i64) :: val
        integer(kind=i64) :: n, i, accumulators(vec_len)
        val = 1_i64
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 1_i64
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators * x(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) * x(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val * accumulators(i) * accumulators(i + 1_i64)
            end do
        end if
    end function prod_i64

    pure function prod_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp) :: val
        real(kind=sp) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 1.0_sp
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 1.0_sp
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators * x(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) * x(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val * accumulators(i) * accumulators(i + 1_i64)
            end do
        end if
    end function prod_sp

    pure function prod_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: val
        real(kind=dp) :: accumulators(vec_len)
        integer(kind=i64) :: n, i
        val = 1.0_dp
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 1.0_dp
            do i = 1_i64, n - vec_len + 1_i64, vec_len
                accumulators = accumulators * x(i:i+vec_len-1_i64)
            end do
            i = n - (n/vec_len)*vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) * x(n-i+1_i64:n)
            do i = 1_i64, vec_len, 2_i64
                val = val * accumulators(i) * accumulators(i + 1_i64)
            end do
        end if
    end function prod_dp

    pure function avg_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp) :: val
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: avg function invalid for array with length < 1')
        val = dsum(x)/real(n, kind=sp)
    end function avg_sp

    pure function avg_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: val
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: avg function invalid for array with length < 1')
        val = dsum(x)/real(n, kind=dp)
    end function avg_dp

    pure function geomean_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp) :: val
        integer(kind=i64) :: n
        call debug_error_condition(logical(any(x <= 0.0_sp), kind=c_bool), &
                                   'module STATISTICS :: geomean function invalid for values <= 0.0')
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: geomean function invalid for array with length < 1')
        val = prod(x)**(1.0_sp/real(n, kind=sp))
    end function geomean_sp

    pure function geomean_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: val
        integer(kind=i64) :: n
        call debug_error_condition(logical(any(x <= 0.0_dp), kind=c_bool), &
                                   'module STATISTICS :: geomean function invalid for values <= 0.0')
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: geomean function invalid for array with length < 1')
        val = prod(x)**(1.0_dp/real(n, kind=dp))
    end function geomean_dp

    pure function std_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp) :: val
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: std function invalid for array with length < 2')
        val = sqrt(dsum((x - avg(x))**2)/real(n - 1_i64, kind=sp))
    end function std_sp

    pure function std_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: val
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: std function invalid for array with length < 2')
        val = sqrt(dsum((x - avg(x))**2)/real(n - 1_i64, kind=dp))
    end function std_dp

    pure function std_weighted_sp(x, w) result(val)
        real(kind=sp), intent(in) :: x(:), w(:)
        real(kind=sp) :: val
        integer(kind=i64) :: n, m
        real(kind=sp) :: w_norm(size(w, kind=i64)), wx_avg, wx_cen2(size(x, kind=i64)), numer, mfrac
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: std function invalid for array with length < 2')
        call debug_error_condition(logical(size(w, kind=i64) /= n, kind=c_bool), &
                                   'module STATISTICS :: weighted std function must have x and w same length')
        call normalize(w, w_norm)
        wx_avg = dsum(w_norm*x)
        wx_cen2 = x - wx_avg !! center x
        wx_cen2 = wx_cen2**2 !! square (x - xbar) term
        numer = dsum(w_norm*wx_cen2) !! numerator under square root
        m = count(.not.nearly(abs(w_norm), 0.0_sp), kind=i64) !! m is quantity of non-zero weights
        mfrac = real(m - 1_i64, kind=sp)/real(m, kind=sp) !! (m-1)/m
        val = sqrt(numer/mfrac)
    end function std_weighted_sp

    pure function std_weighted_dp(x, w) result(val)
        real(kind=dp), intent(in) :: x(:), w(:)
        real(kind=dp) :: val
        integer(kind=i64) :: n, m
        real(kind=dp) :: w_norm(size(w, kind=i64)), wx_avg, wx_cen2(size(x, kind=i64)), numer, mfrac
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: std function invalid for array with length < 2')
        call debug_error_condition(logical(size(w, kind=i64) /= n, kind=c_bool), &
                                   'module STATISTICS :: weighted std function must have x and w same length')
        call normalize(w, w_norm)
        wx_avg = dsum(w_norm*x)
        wx_cen2 = x - wx_avg !! center x
        wx_cen2 = wx_cen2**2 !! square (x - xbar) term
        numer = dsum(w_norm*wx_cen2) !! numerator under square root
        m = count(.not.nearly(abs(w_norm), 0.0_dp), kind=i64) !! m is quantity of non-zero weights
        mfrac = real(m - 1_i64, kind=dp)/real(m, kind=dp) !! (m-1)/m
        val = sqrt(numer/mfrac)
    end function std_weighted_dp

    pure elemental function normal_pdf_1d_sp(x, mu, sig) result(val)
        real(kind=sp), intent(in) :: x, mu, sig
        real(kind=sp) :: val
        real(kind=sp) :: sig2
        call debug_error_condition(sig <= 0.0_sp, &
                                   'module STATISTICS :: normal_pdf function invalid for input with sig == 0.0')
        sig2 = sig*sig
        val = exp(-(x - mu)**2/(2.0_sp*sig2))/sqrt(twopi_sp*sig2)
    end function normal_pdf_1d_sp

    pure elemental function normal_pdf_1d_dp(x, mu, sig) result(val)
        real(kind=dp), intent(in) :: x, mu, sig
        real(kind=dp) :: val
        real(kind=dp) :: sig2
        call debug_error_condition(sig <= 0.0_dp, &
                                   'module STATISTICS :: normal_pdf function invalid for input with sig == 0.0')
        sig2 = sig*sig
        val = exp(-(x - mu)**2/(2.0_dp*sig2))/sqrt(twopi_dp*sig2)
    end function normal_pdf_1d_dp

    pure subroutine normal_pdf_diagonal_covariance_sp(vals, d, x, mu, sig_diag)
        real(kind=sp), intent(out) :: vals(:)
        integer(kind=i64), intent(in) :: d
        real(kind=sp), intent(in) :: x(d,size(vals, kind=i64)), mu(d), sig_diag(d)
        real(kind=sp) :: factor, dx2(d), var_inv(d), var_det
        integer(kind=i64) :: i
        call debug_error_condition(logical(any(sig_diag <= 0.0_sp), kind=c_bool), &
                                   'module STATISTICS :: mv_normal_pdf subroutine invalid for input with det(sig) <= 0.0')
        var_det = prod(sig_diag**2)
        factor = 1.0_sp/sqrt(twopi_sp**d*var_det)
        var_inv = 1.0_sp/sig_diag**2
        do i=1_i64,size(vals, kind=i64)
            dx2 = (x(:,i) - mu)**2
            vals(i) = exp(-0.5_sp*dsum(dx2*var_inv))
        end do
    end subroutine normal_pdf_diagonal_covariance_sp

    pure subroutine normal_pdf_diagonal_covariance_dp(vals, d, x, mu, sig_diag)
        real(kind=dp), intent(out) :: vals(:)
        integer(kind=i64), intent(in) :: d
        real(kind=dp), intent(in) :: x(d,size(vals, kind=i64)), mu(d), sig_diag(d)
        real(kind=dp) :: factor, dx2(d), var_inv(d), var_det
        integer(kind=i64) :: i
        call debug_error_condition(logical(any(sig_diag <= 0.0_dp), kind=c_bool), &
                                   'module STATISTICS :: mv_normal_pdf subroutine invalid for input with det(sig) <= 0.0')
        var_det = prod(sig_diag**2)
        factor = 1.0_dp/sqrt(twopi_dp**d*var_det)
        var_inv = 1.0_dp/sig_diag**2
        do i=1_i64,size(vals, kind=i64)
            dx2 = (x(:,i) - mu)**2
            vals(i) = factor*exp(-0.5_dp*dsum(dx2*var_inv))
        end do
    end subroutine normal_pdf_diagonal_covariance_dp

    pure subroutine normalize_inplace_sp(x)
        real(kind=sp), intent(inout) :: x(:)
        real(kind=sp) :: total
        total = dsum(x)
        call debug_error_condition(nearly(total, 0.0_sp), &
                                   'module STATISTICS :: normalize subroutine invalid for input with total == 0.0')
        x = x/total
    end subroutine normalize_inplace_sp

    pure subroutine normalize_sp(x, xnorm)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp), intent(out) :: xnorm(size(x, kind=i64))
        real(kind=sp) :: total
        total = dsum(x)
        call debug_error_condition(nearly(total, 0.0_sp), &
                                   'module STATISTICS :: normalize subroutine invalid for input with total == 0.0')
        xnorm = x/total
    end subroutine normalize_sp

    pure subroutine normalize_inplace_dp(x)
        real(kind=dp), intent(inout) :: x(:)
        real(kind=dp) :: total
        total = dsum(x)
        call debug_error_condition(nearly(total, 0.0_dp), &
                                   'module STATISTICS :: normalize subroutine invalid for input with total == 0.0')
        x = x/total
    end subroutine normalize_inplace_dp

    pure subroutine normalize_dp(x, xnorm)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(out) :: xnorm(size(x, kind=i64))
        real(kind=dp) :: total
        total = dsum(x)
        call debug_error_condition(nearly(total, 0.0_dp), &
                                   'module STATISTICS :: normalize subroutine invalid for input with total == 0.0')
        xnorm = x/total
    end subroutine normalize_dp

    pure subroutine normalize_min_max_sp(x, xnorm)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp), intent(out) :: xnorm(size(x, kind=i64))
        real(kind=sp) :: xmin, xmax, xdiff
        xmin = minval(x)
        xmax = maxval(x)
        call debug_error_condition(nearly(xmin, xmax), &
                                   'module STATISTICS :: normalize_min_max subroutine invalid for input with min == max')
        xdiff = xmax - xmin
        xnorm = (x - xmin)/xdiff
    end subroutine normalize_min_max_sp

    pure subroutine normalize_min_max_inplace_sp(x)
        real(kind=sp), intent(inout) :: x(:)
        real(kind=sp) :: xmin, xmax, xdiff
        xmin = minval(x)
        xmax = maxval(x)
        call debug_error_condition(nearly(xmin, xmax), &
                                   'module STATISTICS :: normalize_min_max subroutine invalid for input with min == max')
        xdiff = xmax - xmin
        x = (x - xmin)/xdiff
    end subroutine normalize_min_max_inplace_sp

    pure subroutine normalize_min_max_dp(x, xnorm)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(out) :: xnorm(size(x, kind=i64))
        real(kind=dp) :: xmin, xmax, xdiff
        xmin = minval(x)
        xmax = maxval(x)
        call debug_error_condition(nearly(xmin, xmax), &
                                   'module STATISTICS :: normalize_min_max subroutine invalid for input with min == max')
        xdiff = xmax - xmin
        xnorm = (x - xmin)/xdiff
    end subroutine normalize_min_max_dp

    pure subroutine normalize_min_max_inplace_dp(x)
        real(kind=dp), intent(inout) :: x(:)
        real(kind=dp) :: xmin, xmax, xdiff
        xmin = minval(x)
        xmax = maxval(x)
        call debug_error_condition(nearly(xmin, xmax), &
                                   'module STATISTICS :: normalize_min_max subroutine invalid for input with min == max')
        xdiff = xmax - xmin
        x = (x - xmin)/xdiff
    end subroutine normalize_min_max_inplace_dp

    pure subroutine cumsum_i32(x, xcs)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32), intent(out) :: xcs(size(x, kind=i64))
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: cumsum subroutine invalid for array with length < 1')
        xcs(1) = x(1)
        do i=2_i64,n
            xcs(i) = xcs(i-1_i64) + x(i)
        end do
    end subroutine cumsum_i32

    pure subroutine cumsum_inplace_i32(x)
        integer(kind=i32), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: cumsum subroutine invalid for array with length < 1')
        do i=2_i64,n
            x(i) = x(i-1_i64) + x(i)
        end do
    end subroutine cumsum_inplace_i32

    pure subroutine cumsum_i64(x, xcs)
        integer(kind=i64), intent(in) :: x(:)
        integer(kind=i64), intent(out) :: xcs(size(x, kind=i64))
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: cumsum subroutine invalid for array with length < 1')
        xcs(1) = x(1)
        do i=2_i64,n
            xcs(i) = xcs(i-1_i64) + x(i)
        end do
    end subroutine cumsum_i64

    pure subroutine cumsum_inplace_i64(x)
        integer(kind=i64), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: cumsum subroutine invalid for array with length < 1')
        do i=2_i64,n
            x(i) = x(i-1_i64) + x(i)
        end do
    end subroutine cumsum_inplace_i64

    pure subroutine cumsum_sp(x, xcs)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp), intent(out) :: xcs(size(x, kind=i64))
        integer(kind=i64) :: n, i
        real(kind=dp) :: temp
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: cumsum subroutine invalid for array with length < 1')
        temp = real(x(1), kind=dp)
        xcs(1) = real(temp, kind=sp)
        do i=2_i64,n
            temp = temp + real(x(i), kind=dp)
            xcs(i) = real(temp, kind=sp)
        end do
    end subroutine cumsum_sp

    pure subroutine cumsum_inplace_sp(x)
        real(kind=sp), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        real(kind=dp) :: temp
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: cumsum subroutine invalid for array with length < 1')
        temp = real(x(1), kind=dp)
        do i=2_i64,n
            temp = temp + real(x(i), kind=dp)
            x(i) = real(temp, kind=sp)
        end do
    end subroutine cumsum_inplace_sp

    pure subroutine cumsum_dp(x, xcs)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(out) :: xcs(size(x, kind=i64))
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: cumsum subroutine invalid for array with length < 1')
        xcs(1) = x(1)
        do i=2_i64,n
            xcs(i) = xcs(i-1_i64) + x(i)
        end do
    end subroutine cumsum_dp

    pure subroutine cumsum_inplace_dp(x)
        real(kind=dp), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: cumsum subroutine invalid for array with length < 1')
        do i=2_i64,n
            x(i) = x(i-1_i64) + x(i)
        end do
    end subroutine cumsum_inplace_dp

    pure subroutine cov_sp(x, xcov)
        real(kind=sp), intent(in) :: x(:,:)
        real(kind=sp), intent(out) :: xcov(size(x, dim=1, kind=i64),size(x, dim=1, kind=i64))
        integer(kind=i64) :: ndims, i, j
        real(kind=sp) :: x_row(size(x, dim=2, kind=i64)), u, x_centered(size(x, dim=2, kind=i64),size(x, dim=1, kind=i64)), &
                         nsamples_1
        call debug_error_condition(logical(size(x, dim=2, kind=i64) < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: cov subroutine invalid for array with dim=2 < 2')
        ndims = size(x, dim=1, kind=i64)
        do i=1_i64,ndims
            x_row = x(i,:)
            u = avg(x_row)
            x_centered(:,i) = x_row - u
        end do
        nsamples_1 = real(size(x, dim=2, kind=i64) - 1_i64, kind=sp)
        do i=1_i64,ndims
            do j=i,ndims
                xcov(j,i) = vdot(x_centered(:,i), x_centered(:,j))
            end do
            xcov(i:ndims,i) = xcov(i:ndims,i)/nsamples_1
            xcov(i,i) = xcov(i,i) + eps_sp
        end do
        do i=2_i64,ndims
            x_row(1_i64:i-1_i64) = xcov(i,1_i64:i-1_i64)
            do j=1_i64,i-1_i64
                xcov(j,i) = x_row(j)
            end do
        end do
    end subroutine cov_sp

    pure subroutine cov_dp(x, xcov)
        real(kind=dp), intent(in) :: x(:,:)
        real(kind=dp), intent(out) :: xcov(size(x, dim=1, kind=i64),size(x, dim=1, kind=i64))
        integer(kind=i64) :: ndims, i, j
        real(kind=dp) :: x_row(size(x, dim=2, kind=i64)), u, x_centered(size(x, dim=2, kind=i64),size(x, dim=1, kind=i64)), &
                         nsamples_1
        call debug_error_condition(logical(size(x, dim=2, kind=i64) < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: cov subroutine invalid for array with dim=2 < 2')
        ndims = size(x, dim=1, kind=i64)
        do i=1_i64,ndims
            x_row = x(i,:)
            u = avg(x_row)
            x_centered(:,i) = x_row - u
        end do
        nsamples_1 = real(size(x, dim=2, kind=i64) - 1_i64, kind=dp)
        do i=1_i64,ndims
            do j=i,ndims
                xcov(j,i) = vdot(x_centered(:,i), x_centered(:,j))
            end do
            xcov(i:ndims,i) = xcov(i:ndims,i)/nsamples_1
            xcov(i,i) = xcov(i,i) + eps_dp
        end do
        do i=2_i64,ndims
            x_row(1_i64:i-1_i64) = xcov(i,1_i64:i-1_i64)
            do j=1_i64,i-1_i64
                xcov(j,i) = x_row(j)
            end do
        end do
    end subroutine cov_dp

    pure subroutine cov_weighted_sp(x, w, xcov)
        real(kind=sp), intent(in) :: x(:,:), w(:)
        real(kind=sp), intent(out) :: xcov(size(x, dim=1, kind=i64),size(x, dim=1, kind=i64))
        integer(kind=i64) :: ndims, i, j
        real(kind=sp) :: x_row(size(x, dim=2, kind=i64)), u, x_centered(size(x, dim=2, kind=i64),size(x, dim=1, kind=i64)), &
                         bessel, w_norm(size(x, dim=2, kind=i64))
        call debug_error_condition(logical(size(x, dim=2, kind=i64) < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: cov subroutine invalid for array with dim=2 < 2')
        call debug_error_condition(logical(size(x, dim=2, kind=i64) /= size(w), kind=c_bool), &
                                   'module STATISTICS :: cov subroutine number of samples in x must match size of w')
        call normalize(w, w_norm)
        ndims = size(x, dim=1, kind=i64)
        do i=1_i64,ndims
            x_row = x(i,:)
            u = vdot(w_norm, x_row)
            x_centered(:,i) = x_row - u
        end do
        bessel = 1.0_sp/(1.0_sp - vdot(w_norm, w_norm))
        do i=1_i64,ndims
            do j=i,ndims
                xcov(j,i) = dsum(w_norm*x_centered(:,i)*x_centered(:,j))
            end do
            xcov(i:ndims,i) = xcov(i:ndims,i)/bessel
            xcov(i,i) = xcov(i,i) + eps_sp
        end do
        do i=2_i64,ndims
            x_row(1_i64:i-1_i64) = xcov(i,1_i64:i-1_i64)
            do j=1_i64,i-1_i64
                xcov(j,i) = x_row(j)
            end do
        end do
    end subroutine cov_weighted_sp

    pure subroutine cov_weighted_dp(x, w, xcov)
        real(kind=dp), intent(in) :: x(:,:), w(:)
        real(kind=dp), intent(out) :: xcov(size(x, dim=1, kind=i64),size(x, dim=1, kind=i64))
        integer(kind=i64) :: ndims, i, j
        real(kind=dp) :: x_row(size(x, dim=2, kind=i64)), u, x_centered(size(x, dim=2, kind=i64),size(x, dim=1, kind=i64)), &
                         bessel, w_norm(size(x, dim=2, kind=i64))
        call debug_error_condition(logical(size(x, dim=2, kind=i64) < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: cov subroutine invalid for array with dim=2 < 2')
        call debug_error_condition(logical(size(x, dim=2, kind=i64) /= size(w), kind=c_bool), &
                                   'module STATISTICS :: cov subroutine number of samples in x must match size of w')
        call normalize(w, w_norm)
        ndims = size(x, dim=1, kind=i64)
        do i=1_i64,ndims
            x_row = x(i,:)
            u = vdot(w_norm, x_row)
            x_centered(:,i) = x_row - u
        end do
        bessel = 1.0_dp/(1.0_dp - vdot(w_norm, w_norm))
        do i=1_i64,ndims
            do j=i,ndims
                xcov(j,i) = dsum(w_norm*x_centered(:,i)*x_centered(:,j))
            end do
            xcov(i:ndims,i) = xcov(i:ndims,i)/bessel
            xcov(i,i) = xcov(i,i) + eps_dp
        end do
        do i=2_i64,ndims
            x_row(1_i64:i-1_i64) = xcov(i,1_i64:i-1_i64)
            do j=1_i64,i-1_i64
                xcov(j,i) = x_row(j)
            end do
        end do
    end subroutine cov_weighted_dp

end module statistics
