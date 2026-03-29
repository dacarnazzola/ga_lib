module util
use, non_intrinsic :: kinds, only: i32, i64, sp, dp
use, non_intrinsic :: system, only: debug_error_condition
implicit none
private
public :: cov, chol

contains

    pure subroutine cov(c, x, x_avg_opt)
        real(kind=sp), intent(in) :: x(:,:)
        real(kind=sp), intent(out) :: c(size(x, dim=1, kind=i64),size(x, dim=1, kind=i64))
        real(kind=sp), intent(in), optional :: x_avg_opt(:)
        real(kind=dp) :: x_avg(size(x, dim=1, kind=i64)), x_centered(size(x, dim=1, kind=i64),size(x, dim=2, kind=i64)), &
                         c_dp(size(x, dim=1, kind=i64),size(x, dim=1, kind=i64)), inv_n, inv_n_1
        integer(kind=i32) :: nx, i
        call debug_error_condition(size(x,dim=2,kind=i64) > huge(1_i32), &
                                   'UTIL::COV sample x matrix too large for i32 storage')
        call debug_error_condition(size(x,dim=2) < 2_i32, &
                                   'UTIL::COV covariance matrix calculation requires >= 2 samples')
        nx = size(x, dim=2)
        if (present(x_avg_opt)) then
            call debug_error_condition(size(x_avg_opt) /= size(x,dim=1), &
                                       'UTIL::COV average x vector size does not match sample x matrix')
            x_avg = real(x_avg_opt, kind=dp)
        else
            inv_n = 1.0_dp/real(nx, kind=dp)
            x_avg = sum(real(x, kind=dp), dim=2)*inv_n
        end if
        do concurrent (i=1_i32:nx)
            x_centered(:,i) = real(x(:,i), kind=dp) - x_avg
        end do
        inv_n_1 = 1.0_dp/real(nx - 1_i32, kind=dp)
        c_dp = matmul(x_centered, transpose(x_centered))*inv_n_1
        c = real(c_dp, kind=sp)
    end subroutine cov

    pure subroutine chol(l, a)
        real(kind=sp), intent(in) :: a(:,:)
        real(kind=sp), intent(out) :: l(size(a, dim=1, kind=i64),size(a, dim=2, kind=i64))
        real(kind=dp) :: a_dp(size(a, dim=1, kind=i64),size(a, dim=2, kind=i64)), &
                         l_dp(size(a, dim=1, kind=i64),size(a, dim=2, kind=i64)), current_diagonal, sqrt_diag
        integer(kind=i32) :: n, i, j
        call debug_error_condition((size(a,dim=1,kind=i64) > huge(1_i32)) .or. (size(a,dim=2,kind=i64) > huge(1_i32)), &
                                   'UTIL::CHOL supplied a matrix too large for i32 storage')
        call debug_error_condition(size(a,dim=1) /= size(a,dim=2), &
                                   'UTIL::CHOL supplied a matrix must be square')
        n = size(a, dim=1)
        a_dp = real(a, kind=dp)
        l_dp = 0.0_dp
        do i=1_i32,n
            current_diagonal = a_dp(i,i) - dot_product(l_dp(i,1_i32:i-1_i32), l_dp(i,1_i32:i-1_i32))
            if (current_diagonal > 0.0_dp) then
                sqrt_diag = sqrt(current_diagonal)
                l_dp(i,i) = sqrt_diag
                do concurrent (j=i+1_i32:n)
                    l_dp(j,i) = (a_dp(j,i) - dot_product(l_dp(j,1_i32:i-1_i32), l_dp(i,1_i32:i-1_i32)))/sqrt_diag
                end do
            else if (current_diagonal < 0.0_dp) then
                call debug_error_condition(.true., &
                                           'UTIL::CHOL Cholesky decomposition failed, '// &
                                           'supplied a matrix must be at least positive semi-definite')
            end if
        end do
        l = real(l_dp, kind=sp)
    end subroutine chol

end module util
