module crossover
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, bool
use, non_intrinsic :: system, only: debug_error_condition
use, non_intrinsic :: random, only: random_normal_multi
implicit none
private
public :: crossover_fitness_weighted_gaussian_blend

contains

    impure subroutine crossover_fitness_weighted_gaussian_blend(population, selected_pairs_ii, fitness, &
                                                                cholesky_factor, new_population)
        real(kind=sp), intent(in) :: population(:,:), fitness(:), cholesky_factor(:,:)
        integer(kind=i32), intent(in) :: selected_pairs_ii(:,:)
        real(kind=sp), intent(out) :: new_population(:,:)
        real(kind=dp), allocatable :: genes_offspring_spread(:,:), zeros(:), cholesky_factor_dp(:,:)
        integer(kind=i32) :: i
        call debug_error_condition((size(population, dim=1, kind=i64) > huge(1_i32)) .or. &
                                   (size(population, dim=2, kind=i64) > huge(1_i32)) .or. &
                                   (size(selected_pairs_ii, dim=1, kind=i64) > huge(1_i32)) .or. &
                                   (size(selected_pairs_ii, dim=2, kind=i64) > huge(1_i32)) .or. &
                                   (size(fitness, kind=i64) > huge(1_i32)) .or. &
                                   (size(new_population, dim=1, kind=i64) > huge(1_i32)) .or. &
                                   (size(new_population, dim=2, kind=i64) > huge(1_i32)) .or. &
                                   (size(new_population, kind=i64) > huge(1_i32)), &
                                   'CROSSOVER::crossover_fitness_weighted_gaussian_blend input arrays too large for i32 storage')
        call debug_error_condition(size(population, dim=2) /= size(fitness), &
                                   'CROSSOVER::crossover_fitness_weighted_gaussian_blend population size must match fitness size')
        call debug_error_condition(size(selected_pairs_ii, dim=1) /= 2, &
                                   'CROSSOVER::crossover_fitness_weighted_gaussian_blend crossover only implemented for 2 parents')
        call debug_error_condition(size(selected_pairs_ii, dim=2) /= size(new_population, dim=2), &
                                   'CROSSOVER::crossover_fitness_weighted_gaussian_blend crossover pairs mismatch new population')
        call debug_error_condition(size(population, dim=1) /= size(new_population, dim=1), &
                                   'CROSSOVER::crossover_fitness_weighted_gaussian_blend new population genes mismatch original')
        call debug_error_condition((minval(selected_pairs_ii) < 1) .or. &
                                   (maxval(selected_pairs_ii) > size(population, dim=2)), &
                                   'CROSSOVER::crossover_fitness_weighted_gaussian_blend crossover pairs exceed population indices')
        allocate(genes_offspring_spread(size(new_population, dim=1),size(new_population, dim=2)), &
                 zeros(size(new_population, dim=1)), cholesky_factor_dp(size(new_population, dim=1),size(new_population, dim=1)))
        zeros = 0.0_dp
        cholesky_factor_dp = real(cholesky_factor, kind=dp)
        call random_normal_multi(genes_offspring_spread, zeros, cholesky_factor_dp)
        do concurrent (i=1_i32:size(new_population, dim=2))
            block
                real(kind=dp) :: genes1(size(population, dim=1)), genes2(size(population, dim=1)), &
                                 fitness1, fitness2, sum_fitness, w1, w2, genes_scale
                logical(kind=bool) :: good_sum

                genes1 = real(population(:,selected_pairs_ii(1,i)), kind=dp)
                genes2 = real(population(:,selected_pairs_ii(2,i)), kind=dp)
                fitness1 = real(fitness(selected_pairs_ii(1,i)), kind=dp)
                fitness2 = real(fitness(selected_pairs_ii(2,i)), kind=dp)
                sum_fitness = fitness1 + fitness2
                good_sum = sum_fitness >= 1.0e-30_dp
                w1 = merge(fitness2/sum_fitness, 0.5_dp, good_sum)
                w2 = merge(fitness1/sum_fitness, 0.5_dp, good_sum)
                genes_scale = min(norm2(w1*genes1 + w2*genes2 - genes1), &
                                  norm2(w1*genes1 + w2*genes2 - genes2))/real(size(population, dim=1), kind=dp)

                new_population(:,i) = real(w1*genes1 + w2*genes2 + genes_scale*genes_offspring_spread(:,i), kind=sp)
            end block
        end do
    end subroutine crossover_fitness_weighted_gaussian_blend

end module crossover
