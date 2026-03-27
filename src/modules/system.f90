module system
use, non_intrinsic :: kinds, only: sp, dp, c_bool, c_char
use, non_intrinsic :: constants, only: debug
implicit none
private

    interface debug_error_condition
        module procedure :: debug_error_condition_default
        module procedure :: debug_error_condition_c_bool
    end interface debug_error_condition

    interface nearly
        module procedure :: nearly_sp
        module procedure :: nearly_dp
    end interface nearly

    public :: debug_error_condition, nearly

contains

    pure elemental subroutine debug_error_condition_default(logical_expression, error_message)
        logical, intent(in) :: logical_expression
        character(len=*, kind=c_char), intent(in) :: error_message
        if (debug) then
            if (logical_expression) then
                error stop error_message
            end if
        end if
    end subroutine debug_error_condition_default

    pure elemental subroutine debug_error_condition_c_bool(logical_expression, error_message)
        logical(kind=c_bool), intent(in) :: logical_expression
        character(len=*, kind=c_char), intent(in) :: error_message
        if (debug) then
            if (logical_expression) then
                error stop error_message
            end if
        end if
    end subroutine debug_error_condition_c_bool

    pure elemental function nearly_sp(x, check) result(val)
        real(kind=sp), intent(in) :: x, check
        logical(kind=c_bool) :: val
        val = logical(abs(check - x) < spacing(check), kind=c_bool)
    end function nearly_sp

    pure elemental function nearly_dp(x, check) result(val)
        real(kind=dp), intent(in) :: x, check
        logical(kind=c_bool) :: val
        val = logical(abs(check - x) < spacing(check), kind=c_bool)
    end function nearly_dp

end module system
