module mutation
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, bool
use, non_intrinsic :: system, only: debug_error_condition
use, non_intrinsic :: random, only: random_uniform_dp, random_normal_multi
implicit none
private
public :: mutation_covariance_aware_gaussian

contains

    impure subroutine mutation_covariance_aware_gaussian(population, mutation_rate, cholesky_factor)
        real(kind=sp), intent(inout) :: population(:,:)
        real(kind=sp), intent(in) :: mutation_rate, cholesky_factor(:,:)
        real(kind=dp), allocatable :: do_mutate(:), mutation_magnitude(:,:), cholesky_factor_dp(:,:), zeros(:)
        integer(kind=i32), allocatable :: mutate_ii(:)
        logical(kind=bool), allocatable :: pass_mutate(:)
        integer(kind=i32) :: i, num_mutate
        call debug_error_condition((size(population, dim=1, kind=i64) > huge(1_i32)) .or. &
                                   (size(population, dim=2, kind=i64) > huge(1_i32)) .or. &
                                   (size(cholesky_factor, dim=1, kind=i64) > huge(1_i32)) .or. &
                                   (size(cholesky_factor, dim=2, kind=i64) > huge(1_i32)), &
                                   'MUTATION::MUTATION_COVARIANCE_AWARE_GAUSSIAN input matrices too large for i32 storage')
        call debug_error_condition((mutation_rate < 0.0_dp) .or. (mutation_rate > 1.0_dp), &
                                   'MUTATION::MUTATION_COVARIANCE_AWARE_GAUSSIAN mutation rate must be [0.0, 1.0]')
        call debug_error_condition(size(population, dim=1) /= size(cholesky_factor, dim=1), &
                                   'MUTATION::MUTATION_COVARIANCE_AWARE_GAUSSIAN Cholesky factor must match gene dimension')
        call debug_error_condition(size(cholesky_factor, dim=1) /= size(cholesky_factor, dim=2), &
                                   'MUTATION::MUTATION_COVARIANCE_AWARE_GAUSSIAN Cholesky factor must be a square matrix')
        call debug_error_condition(any([(cholesky_factor(i,i), i=1,size(cholesky_factor, dim=1))] < 0.0_sp), &
                                   'MUTATION::MUTATION_COVARIANCE_AWARE_GAUSSIAN Cholesky factor malformed with negative diagonal')
        allocate(do_mutate(size(population, dim=2)), pass_mutate(size(population, dim=2)))
        call random_uniform_dp(do_mutate, size(do_mutate, kind=i32), 0.0_dp, 1.0_dp)
        pass_mutate = logical(do_mutate < real(mutation_rate, kind=dp), kind=bool)
        num_mutate = count(pass_mutate, kind=i32)
        if (num_mutate > 0_i32) then
            allocate(mutate_ii(num_mutate), &
                     mutation_magnitude(size(population, dim=1),num_mutate), &
                     cholesky_factor_dp(size(cholesky_factor, dim=1),size(cholesky_factor, dim=2)), &
                     zeros(size(population, dim=1)))
            cholesky_factor_dp = real(cholesky_factor, kind=dp)
            zeros = 0.0_dp
            mutate_ii = pack([(i, i=1,size(pass_mutate))], pass_mutate)
            call random_normal_multi(mutation_magnitude, zeros, cholesky_factor_dp)
            do concurrent (i=1_i32:num_mutate)
                population(:,mutate_ii(i)) = real(real(population(:,mutate_ii(i)), kind=dp) + mutation_magnitude(:,i), kind=sp)
            end do
        end if
    end subroutine mutation_covariance_aware_gaussian

end module mutation
