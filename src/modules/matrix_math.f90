module matrix_math
use, non_intrinsic :: kinds, only: i64, sp, dp, c_bool
use, non_intrinsic :: constants, only: eps_sp, eps_dp
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: vector_math, only: vmag, vdot
implicit none
private

    interface chol
        module procedure :: chol_sp
        module procedure :: chol_dp
    end interface chol

    interface qr
        module procedure :: qr_opt_tau_sp
        module procedure :: qr_opt_tau_dp
    end interface qr

    interface forward_substitution
        module procedure :: forward_substitution_sp
        module procedure :: forward_substitution_dp
    end interface forward_substitution

    interface backward_substitution
        module procedure :: backward_substitution_sp
        module procedure :: backward_substitution_dp
    end interface backward_substitution

    public :: chol, qr, forward_substitution, backward_substitution

contains

    pure subroutine chol_sp(A, L)
        real(kind=sp), intent(in) :: A(:,:)
        real(kind=sp), intent(out) :: L(size(A, dim=1, kind=i64),size(A, dim=2, kind=i64))
        integer(kind=i64) :: n, i, j, k
        real(kind=sp) :: sum_of_squares, tol
        call debug_error_condition(logical(size(A, dim=1, kind=i64) /= size(A, dim=2, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine chol requires A is a square (n x n) matrix')
        tol = max(eps_sp, eps_sp*maxval(abs(A)))
        L = 0.0_sp
        n = size(A, dim=1, kind=i64)
        do i=1_i64,n
            sum_of_squares = A(i,i)
            do j=1_i64,i-1_i64
                sum_of_squares = sum_of_squares - L(i,j)**2
            end do
            if (sum_of_squares > tol) then
                L(i,i) = sqrt(sum_of_squares)
                do j=i+1_i64,n
                    sum_of_squares = A(j,i)
                    do k=1_i64,i-1_i64
                        sum_of_squares = sum_of_squares - L(j,k)*L(i,k)
                    end do
                    L(j,i) = sum_of_squares/L(i,i)
                end do
            else if (sum_of_squares > -tol) then
                L(:,i) = 0.0_sp
            else
                call debug_error_condition(logical(.true., kind=c_bool), &
                                           'module MATRIX_MATH :: subroutine chol requires A be positive semi-definite')
            end if
        end do
    end subroutine chol_sp

    pure subroutine chol_dp(A, L)
        real(kind=dp), intent(in) :: A(:,:)
        real(kind=dp), intent(out) :: L(size(A, dim=1, kind=i64),size(A, dim=2, kind=i64))
        integer(kind=i64) :: n, i, j, k
        real(kind=dp) :: sum_of_squares, tol
        call debug_error_condition(logical(size(A, dim=1, kind=i64) /= size(A, dim=2, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine chol requires A is a square (n x n) matrix')
        tol = max(eps_dp, eps_dp*maxval(abs(A)))
        L = 0.0_dp
        n = size(A, dim=1, kind=i64)
        do i=1_i64,n
            sum_of_squares = A(i,i)
            do j=1_i64,i-1_i64
                sum_of_squares = sum_of_squares - L(i,j)**2
            end do
            if (sum_of_squares > tol) then
                L(i,i) = sqrt(sum_of_squares)
                do j=i+1_i64,n
                    sum_of_squares = A(j,i)
                    do k=1_i64,i-1_i64
                        sum_of_squares = sum_of_squares - L(j,k)*L(i,k)
                    end do
                    L(j,i) = sum_of_squares/L(i,i)
                end do
            else if (sum_of_squares > -tol) then
                L(:,i) = 0.0_dp
            else
                call debug_error_condition(logical(.true., kind=c_bool), &
                                           'module MATRIX_MATH :: subroutine chol requires A be positive semi-definite')
            end if
        end do
    end subroutine chol_dp

    pure subroutine qr_opt_tau_sp(A, opt_tau)
        real(kind=sp), intent(inout) :: A(:,:)
        real(kind=sp), intent(out), optional :: opt_tau(min(size(A, dim=1, kind=i64), size(A, dim=2, kind=i64)))
        real(kind=sp) :: x(size(A, dim=1, kind=i64)), s, alpha, p(size(A, dim=1, kind=i64)), unorm2, &
                         tau(min(size(A, dim=1, kind=i64), size(A, dim=2, kind=i64)))
        integer(kind=i64) :: m, n, imax, i, x_dim, p_dim
        m = size(A, dim=1, kind=i64)
        n = size(A, dim=2, kind=i64)
        imax = min(m - 1_i64, n)
        do i=1,imax
            x_dim = m - i + 1_i64
            x(1_i64:x_dim) = A(i:m,i)
            s = vmag(x(1_i64:x_dim))
            if (nearly(s, 0.0_sp)) then
                tau(i) = 0.0_sp
                cycle
            end if
            alpha = -1.0_sp*sign(s, x(1_i64))
            x(1_i64) = x(1_i64) - alpha
            unorm2 = vdot(x(1_i64:x_dim), x(1_i64:x_dim))
            if (nearly(unorm2, 0.0_sp)) then
                tau(i) = 0.0_sp
                cycle
            end if
            tau(i) = 2.0_sp*x(1_i64)**2/unorm2
            x(1_i64:x_dim) = x(1_i64:x_dim)/x(1_i64)
            p_dim = n - i + 1_i64
            p(1_i64:p_dim) = reshape(matmul(reshape(x(1_i64:x_dim), shape=[1_i64, x_dim]), A(i:m,i:n)), shape=[p_dim])
            A(i:m,i:n) = A(i:m,i:n) - tau(i)*spread(x(1_i64:x_dim), dim=2, ncopies=p_dim) * &
                                             spread(p(1_i64:p_dim), dim=1, ncopies=x_dim)
            A(i+1_i64:m,i) = x(2_i64:m-i+1_i64)
        end do
        if (m <= n) tau(m) = 0.0_sp
        if (present(opt_tau)) opt_tau = tau
    end subroutine qr_opt_tau_sp

    pure subroutine qr_opt_tau_dp(A, opt_tau)
        real(kind=dp), intent(inout) :: A(:,:)
        real(kind=dp), intent(out), optional :: opt_tau(min(size(A, dim=1, kind=i64), size(A, dim=2, kind=i64)))
        real(kind=dp) :: x(size(A, dim=1, kind=i64)), s, alpha, p(size(A, dim=1, kind=i64)), unorm2, &
                         tau(min(size(A, dim=1, kind=i64), size(A, dim=2, kind=i64)))
        integer(kind=i64) :: m, n, imax, i, x_dim, p_dim
        m = size(A, dim=1, kind=i64)
        n = size(A, dim=2, kind=i64)
        imax = min(m - 1_i64, n)
        do i=1,imax
            x_dim = m - i + 1_i64
            x(1_i64:x_dim) = A(i:m,i)
            s = vmag(x(1_i64:x_dim))
            if (nearly(s, 0.0_dp)) then
                tau(i) = 0.0_dp
                cycle
            end if
            alpha = -1.0_dp*sign(s, x(1_i64))
            x(1_i64) = x(1_i64) - alpha
            unorm2 = vdot(x(1_i64:x_dim), x(1_i64:x_dim))
            if (nearly(unorm2, 0.0_dp)) then
                tau(i) = 0.0_dp
                cycle
            end if
            tau(i) = 2.0_dp*x(1_i64)**2/unorm2
            x(1_i64:x_dim) = x(1_i64:x_dim)/x(1_i64)
            p_dim = n - i + 1_i64
            p(1_i64:p_dim) = reshape(matmul(reshape(x(1_i64:x_dim), shape=[1_i64, x_dim]), A(i:m,i:n)), shape=[p_dim])
            A(i:m,i:n) = A(i:m,i:n) - tau(i)*spread(x(1_i64:x_dim), dim=2, ncopies=p_dim) * &
                                             spread(p(1_i64:p_dim), dim=1, ncopies=x_dim)
            A(i+1_i64:m,i) = x(2_i64:m-i+1_i64)
        end do
        if (m <= n) tau(m) = 0.0_dp
        if (present(opt_tau)) opt_tau = tau
    end subroutine qr_opt_tau_dp

    pure subroutine forward_substitution_sp(L, b, x)
        real(kind=sp), intent(in) :: L(:,:), b(:)
        real(kind=sp), intent(out) :: x(:)
        integer(kind=i64) :: i
        real(kind=sp) :: temp
        call debug_error_condition(logical(size(b, kind=i64) < 1_i64, kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires size b vector >= 1')
        call debug_error_condition(logical(size(L, dim=1, kind=i64) /= size(L, dim=2, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires L be square (n x n) matrix')
        call debug_error_condition(logical(any([(nearly(L(i,i), 0.0_sp), i=1_i64,size(L,dim=1,kind=i64))]), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires L(i,i) NOT nearly 0.0')
        call debug_error_condition(logical(size(b, kind=i64) /= size(L, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires b vector same dimension as L')
        call debug_error_condition(logical(size(x, kind=i64) /= size(L, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires x vector same dimension as L')
        x(1_i64) = b(1_i64)/L(1_i64,1_i64)
        do i=2_i64,size(b, kind=i64)
            temp = vdot(L(i,1_i64:i-1_i64), x(1_i64:i-1_i64))
            x(i) = (b(i) - temp)/L(i,i)
        end do
    end subroutine forward_substitution_sp

    pure subroutine forward_substitution_dp(L, b, x)
        real(kind=dp), intent(in) :: L(:,:), b(:)
        real(kind=dp), intent(out) :: x(:)
        integer(kind=i64) :: i
        real(kind=dp) :: temp
        call debug_error_condition(logical(size(b, kind=i64) < 1_i64, kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires size b vector >= 1')
        call debug_error_condition(logical(size(L, dim=1, kind=i64) /= size(L, dim=2, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires L be square (n x n) matrix')
        call debug_error_condition(logical(any([(nearly(L(i,i), 0.0_dp), i=1_i64,size(L,dim=1,kind=i64))]), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires L(i,i) NOT nearly 0.0')
        call debug_error_condition(logical(size(b, kind=i64) /= size(L, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires b vector same dimension as L')
        call debug_error_condition(logical(size(x, kind=i64) /= size(L, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine forward_substitution requires x vector same dimension as L')
        x(1_i64) = b(1_i64)/L(1_i64,1_i64)
        do i=2_i64,size(b, kind=i64)
            temp = vdot(L(i,1_i64:i-1_i64), x(1_i64:i-1_i64))
            x(i) = (b(i) - temp)/L(i,i)
        end do
    end subroutine forward_substitution_dp

    pure subroutine backward_substitution_sp(U, b, x)
        real(kind=sp), intent(in) :: U(:,:), b(:)
        real(kind=sp), intent(out) :: x(:)
        integer(kind=i64) :: i, n
        real(kind=sp) :: temp
        call debug_error_condition(logical(size(b, kind=i64) < 1_i64, kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires size b vector >= 1')
        call debug_error_condition(logical(size(U, dim=1, kind=i64) /= size(U, dim=2, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires U be square (n x n) matrix')
        call debug_error_condition(logical(any([(nearly(U(i,i), 0.0_sp), i=1_i64,size(U,dim=1,kind=i64))]), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires U(i,i) NOT nearly 0.0')
        call debug_error_condition(logical(size(b, kind=i64) /= size(U, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires b vector same dimension as U')
        call debug_error_condition(logical(size(x, kind=i64) /= size(U, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires x vector same dimension as U')
        n = size(b, kind=i64)
        x(n) = b(n)/U(n,n)
        do i=n-1_i64,1_i64,-1_i64
            temp = vdot(U(i,i+1_i64:n), x(i+1_i64:n))
            x(i) = (b(i) - temp)/U(i,i)
        end do
    end subroutine backward_substitution_sp

    pure subroutine backward_substitution_dp(U, b, x)
        real(kind=dp), intent(in) :: U(:,:), b(:)
        real(kind=dp), intent(out) :: x(:)
        integer(kind=i64) :: i, n
        real(kind=dp) :: temp
        call debug_error_condition(logical(size(b, kind=i64) < 1_i64, kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires size b vector >= 1')
        call debug_error_condition(logical(size(U, dim=1, kind=i64) /= size(U, dim=2, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires U be square (n x n) matrix')
        call debug_error_condition(logical(any([(nearly(U(i,i), 0.0_dp), i=1_i64,size(U,dim=1,kind=i64))]), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires U(i,i) NOT nearly 0.0')
        call debug_error_condition(logical(size(b, kind=i64) /= size(U, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires b vector same dimension as U')
        call debug_error_condition(logical(size(x, kind=i64) /= size(U, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine backward_substitution requires x vector same dimension as U')
        n = size(b, kind=i64)
        x(n) = b(n)/U(n,n)
        do i=n-1_i64,1_i64,-1_i64
            temp = vdot(U(i,i+1_i64:n), x(i+1_i64:n))
            x(i) = (b(i) - temp)/U(i,i)
        end do
    end subroutine backward_substitution_dp

end module matrix_math
