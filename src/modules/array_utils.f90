module array_utils
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, c_bool
implicit none
private

    interface swap
        module procedure :: swap_i32
        module procedure :: swap_i64
        module procedure :: swap_sp
        module procedure :: swap_dp
        module procedure :: swap_c_bool
    end interface swap

    interface flip
        module procedure :: flip_i32
        module procedure :: flip_i64
        module procedure :: flip_sp
        module procedure :: flip_dp
        module procedure :: flip_c_bool
    end interface flip

    interface broadcast_add
        module procedure :: broadcast_add_i32
        module procedure :: broadcast_add_i64
        module procedure :: broadcast_add_sp
        module procedure :: broadcast_add_dp
    end interface broadcast_add

    interface broadcast_sub
        module procedure :: broadcast_sub_i32
        module procedure :: broadcast_sub_i64
        module procedure :: broadcast_sub_sp
        module procedure :: broadcast_sub_dp
        module procedure :: broadcast_sub_inplace_dp
    end interface broadcast_sub

    interface broadcast_mul
        module procedure :: broadcast_mul_i32
        module procedure :: broadcast_mul_i64
        module procedure :: broadcast_mul_sp
        module procedure :: broadcast_mul_dp
    end interface broadcast_mul

    interface broadcast_div
        module procedure :: broadcast_div_i32
        module procedure :: broadcast_div_i64
        module procedure :: broadcast_div_sp
        module procedure :: broadcast_div_dp
    end interface broadcast_div

    public :: swap, flip, broadcast_add, broadcast_sub, broadcast_mul, broadcast_div

contains

    pure elemental subroutine swap_i32(a, b)
        integer(kind=i32), intent(inout) :: a, b
        integer(kind=i32) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_i32

    pure elemental subroutine swap_i64(a, b)
        integer(kind=i64), intent(inout) :: a, b
        integer(kind=i64) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_i64

    pure elemental subroutine swap_sp(a, b)
        real(kind=sp), intent(inout) :: a, b
        real(kind=sp) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_sp

    pure elemental subroutine swap_dp(a, b)
        real(kind=dp), intent(inout) :: a, b
        real(kind=dp) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_dp

    pure elemental subroutine swap_c_bool(a, b)
        logical(kind=c_bool), intent(inout) :: a, b
        logical(kind=c_bool) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_c_bool

    pure subroutine flip_i32(x)
        integer(kind=i32), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_i32

    pure subroutine flip_i64(x)
        integer(kind=i64), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_i64

    pure subroutine flip_sp(x)
        real(kind=sp), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_sp

    pure subroutine flip_dp(x)
        real(kind=dp), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_dp

    pure subroutine flip_c_bool(x)
        logical(kind=c_bool), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_c_bool

    pure subroutine broadcast_add_i32(arr, x_spread, vals)
        integer(kind=i32), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i32), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i) + x_spread
        end do
    end subroutine broadcast_add_i32

    pure subroutine broadcast_add_i64(arr, x_spread, vals)
        integer(kind=i64), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i64), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i) + x_spread
        end do
    end subroutine broadcast_add_i64

    pure subroutine broadcast_add_sp(arr, x_spread, vals)
        real(kind=sp), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        real(kind=sp), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i) + x_spread
        end do
    end subroutine broadcast_add_sp

    pure subroutine broadcast_add_dp(arr, x_spread, vals)
        real(kind=dp), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        real(kind=dp), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i) + x_spread
        end do
    end subroutine broadcast_add_dp

    pure subroutine broadcast_sub_i32(arr, x_spread, vals)
        integer(kind=i32), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i32), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i) - x_spread
        end do
    end subroutine broadcast_sub_i32

    pure subroutine broadcast_sub_i64(arr, x_spread, vals)
        integer(kind=i64), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i64), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i) - x_spread
        end do
    end subroutine broadcast_sub_i64

    pure subroutine broadcast_sub_sp(arr, x_spread, vals)
        real(kind=sp), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        real(kind=sp), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i) - x_spread
        end do
    end subroutine broadcast_sub_sp

    pure subroutine broadcast_sub_dp(arr, x_spread, vals)
        real(kind=dp), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        real(kind=dp), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i) - x_spread
        end do
    end subroutine broadcast_sub_dp

    pure subroutine broadcast_sub_inplace_dp(arr, x_spread)
        real(kind=dp), intent(inout) :: arr(:,:)
        real(kind=dp), intent(in) :: x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            arr(:,i) = arr(:,i) - x_spread
        end do
    end subroutine broadcast_sub_inplace_dp

    pure subroutine broadcast_mul_i32(arr, x_spread, vals)
        integer(kind=i32), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i32), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i)*x_spread
        end do
    end subroutine broadcast_mul_i32

    pure subroutine broadcast_mul_i64(arr, x_spread, vals)
        integer(kind=i64), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i64), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i)*x_spread
        end do
    end subroutine broadcast_mul_i64

    pure subroutine broadcast_mul_sp(arr, x_spread, vals)
        real(kind=sp), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        real(kind=sp), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i)*x_spread
        end do
    end subroutine broadcast_mul_sp

    pure subroutine broadcast_mul_dp(arr, x_spread, vals)
        real(kind=dp), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        real(kind=dp), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i)*x_spread
        end do
    end subroutine broadcast_mul_dp

    pure subroutine broadcast_div_i32(arr, x_spread, vals)
        integer(kind=i32), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i32), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i)/x_spread
        end do
    end subroutine broadcast_div_i32

    pure subroutine broadcast_div_i64(arr, x_spread, vals)
        integer(kind=i64), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        integer(kind=i64), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i)/x_spread
        end do
    end subroutine broadcast_div_i64

    pure subroutine broadcast_div_sp(arr, x_spread, vals)
        real(kind=sp), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        real(kind=sp), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i)/x_spread
        end do
    end subroutine broadcast_div_sp

    pure subroutine broadcast_div_dp(arr, x_spread, vals)
        real(kind=dp), intent(in) :: arr(:,:), x_spread(size(arr, dim=1, kind=i64))
        real(kind=dp), intent(out) :: vals(size(arr, dim=1, kind=i64),size(arr, dim=2, kind=i64))
        integer(kind=i64) :: i
        do concurrent (i=1_i64:size(arr, dim=2, kind=i64))
            vals(:,i) = arr(:,i)/x_spread
        end do
    end subroutine broadcast_div_dp

end module array_utils
