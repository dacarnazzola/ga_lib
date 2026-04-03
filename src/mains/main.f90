module ga_interface
use, non_intrinsic :: kinds, only: stdout, ik=>i32, rk=>sp
use, non_intrinsic :: constants, only: twopi=>twopi_sp
use, non_intrinsic :: util, only: covariance=>cov, cholesky_decomposition=>chol, sort_candidates=>sort, &
                                  apply_constraints=>constraints_reflective_boundary
use, non_intrinsic :: random, only: random_uniform=>random_uniform_sp
use, non_intrinsic :: selection, only: perform_selection=>selection_tournament
use, non_intrinsic :: crossover, only: perform_crossover=>crossover_fitness_weighted_gaussian_blend
use, non_intrinsic :: mutation, only: perform_mutation=>mutation_covariance_aware_gaussian
implicit none
private
public :: stdout, ik, rk
public :: twopi
public :: covariance, cholesky_decomposition, sort_candidates, apply_constraints, random_uniform
public :: perform_selection, perform_crossover, perform_mutation
end module ga_interface                                                           

module benchmark
use, non_intrinsic :: ga_interface
implicit none
private
public :: ik, solve_rastrigin

contains
    
    pure subroutine rastrigin(x, fx)
        real(rk), intent(in) :: x(:,:)
        real(rk), intent(out) :: fx(:)
        integer :: i, xdim, nx
        xdim = size(x, dim=1)
        nx = size(x, dim=2)
        do concurrent (i=1:nx)
            fx(i) = 10.0_rk*xdim + sum(x(:,i)*x(:,i) - 10.0_rk*cos(twopi*x(:,i)), dim=1)
        end do
    end subroutine rastrigin

    impure subroutine solve_rastrigin(problem_dimension, population_size, maximum_generations)
        integer(ik), intent(in) :: problem_dimension, population_size, maximum_generations
        real(rk), allocatable :: fitness(:), cov(:,:), chol(:,:), baseline(:,:), baseline_fitness(:), domain_lb(:), domain_ub(:), &
                                 regularization_vector(:), candidates(:,:), candidate_fitness(:)
        real(rk), allocatable, target :: pop1(:,:), pop2(:,:)
        real(rk), pointer :: current_population(:,:), new_population(:,:), dummy_ptr(:,:)
        real(rk) :: mutation_rate, mutation_scale, ft, ft_new
        integer(ik), allocatable :: selected_pairs_ii(:,:), candidate_sorted_ii(:)
        integer(ik) :: generation, elite_ii, failed_gen, i, total_evals

        !allocate arrays
        allocate(fitness(population_size), cov(problem_dimension,problem_dimension), chol(problem_dimension,problem_dimension), &
                 domain_lb(problem_dimension), domain_ub(problem_dimension), regularization_vector(problem_dimension), &
                 candidates(problem_dimension,2*population_size), candidate_fitness(2*population_size), &
                 pop1(problem_dimension,population_size), pop2(problem_dimension,population_size), &
                 selected_pairs_ii(2,population_size), candidate_sorted_ii(2*population_size))

        ! initialize population
        domain_lb = -5.12_rk
        domain_ub = 5.12_rk
        call random_uniform(pop1, size(pop1), minval(domain_lb), maxval(domain_ub))
        current_population => pop1
        new_population => pop2

        ! calculate fitness
        call rastrigin(current_population, fitness)
        total_evals = size(fitness)
        elite_ii = minloc(fitness, dim=1)
        write(stdout,'(a,f0.6)') 'initial best fitness: ',fitness(elite_ii)
        ft = fitness(1)

        ! set regularization vector very small, just to avoid numerical collapse
        regularization_vector = (1.0e-4_rk)**2
        ! set mutation rate as 1.0 - 4.0/population_size, enabling high mutation rate for populations 10+
        mutation_rate = 1.0_rk - 4.0_rk/real(population_size, kind=rk)
        ! start mutation scale at 1.0, it will vary depending on generational fitness
        mutation_scale = 1.0_rk

        failed_gen = 0
        ! perform GA loop
        do generation=1,maximum_generations
            ! tournament selection, K=2
            call perform_selection(current_population, fitness, selected_pairs_ii, k_opt=2)

            ! store current generation fitness scores
            candidate_fitness(1:population_size) = fitness

            ! fitness-weighted crossover with Gaussian blend
            call perform_crossover(current_population, selected_pairs_ii, fitness, new_population)
            call apply_constraints(new_population, domain_lb, domain_ub)

            ! calculate covariance matrix and Cholesky factor for mutation
            call covariance(cov, current_population(:,1:population_size/4), reg_vec_opt=regularization_vector)
            call cholesky_decomposition(chol, cov)

            ! Gaussian mutation based on post-crossover population genetic covariance
            call perform_mutation(new_population, mutation_rate, chol, mutation_scale)
            call apply_constraints(new_population, domain_lb, domain_ub)

            ! calculate fitness
            call rastrigin(new_population, fitness)
            total_evals = total_evals + size(fitness)

            ! store current_population and new_population into candidate_population(2*population_size), then keep top half
            candidates(:,1:population_size) = current_population
            candidates(:,population_size+1:2*population_size) = new_population
            candidate_fitness(population_size+1:2*population_size) = fitness
            call sort_candidates(candidate_fitness, candidate_sorted_ii)
            fitness = candidate_fitness(1:population_size)
            new_population = candidates(:,candidate_sorted_ii(1:population_size))

            elite_ii = 1
            write(stdout,'(a,i0,a,f0.6)') 'generation: ',generation,', best fitness: ',fitness(elite_ii)
            ft_new = fitness(1)

            if (ft_new < ft) then
                mutation_scale = max(0.1_rk*mutation_scale, 1.0e-5_rk)
                failed_gen = 0
            else
                mutation_scale = 2.0_rk*mutation_scale
                failed_gen = failed_gen + 1
            end if

            if (failed_gen < 10) then
                ft = ft_new
            else ! elite 1 failing to improve, reset search space
                chol = 0.0_rk
                do concurrent (i=1:problem_dimension)
                    chol(i,i) = abs(domain_ub(i) - domain_lb(i))/2.0_rk
                end do
                do concurrent (i=population_size/2:population_size)
                    new_population(:,i) = new_population(:,1)
                end do
                call perform_mutation(new_population(:,population_size/2:population_size), 1.0_rk, chol, 1.0_rk)
                call apply_constraints(new_population(:,population_size/2:population_size), domain_lb, domain_ub)
                call rastrigin(new_population(:,population_size/2:population_size), fitness(population_size/2:population_size))
                total_evals = total_evals + (population_size - population_size/2)
                elite_ii = minloc(fitness, dim=1)
                write(stdout,'(a,i0,a,f0.6)') 'CATASTROPHE generation: ',generation,', best fitness: ',fitness(elite_ii)
                ft = sum(fitness)
                mutation_scale = 1.0_rk ! reset mutation_scale to 1.0 for randomized population
                failed_gen = 0
            end if

            ! swap population pointers
            dummy_ptr => current_population
            current_population => new_population
            new_population => dummy_ptr
        end do

        ! establish baseline
        allocate(baseline(problem_dimension,total_evals), baseline_fitness(total_evals))
        call random_uniform(baseline, size(baseline), -5.12, 5.12)
        call rastrigin(baseline, baseline_fitness)
        elite_ii = minloc(baseline_fitness, dim=1)
        write(stdout,'(a,f0.6,a,i0,a)') 'baseline best fitness: ',baseline_fitness(elite_ii),' (',total_evals,' evaluations)'
    end subroutine

end module benchmark

program main
use benchmark
implicit none

    integer(ik), parameter :: problem_dimension   = 20
    integer(ik), parameter :: population_size     = 128*9
    integer(ik), parameter :: maximum_generations = 24*3

    call solve_rastrigin(problem_dimension, population_size, maximum_generations)

end program main
