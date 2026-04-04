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
public :: rastrigin, rosenbrock
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
    pure subroutine rosenbrock(x, fx)
        real(rk), intent(in) :: x(:,:)
        real(rk), intent(out) :: fx(:)
        integer :: i, xdim, nx
        xdim = size(x, dim=1)
        nx = size(x, dim=2)
        do concurrent (i=1:nx)
            fx(i) = sum(100.0_rk*(x(2:xdim,i) - x(1:xdim-1,i)**2)**2 + (1.0_rk - x(1:xdim-1,i))**2, dim=1)
        end do
    end subroutine rosenbrock
end module ga_interface                                                           

module benchmark
use, non_intrinsic :: ga_interface, evaluate_function=>rastrigin
implicit none
private
public :: ik, solve

contains

    impure subroutine solve(problem_dimension, population_size, maximum_generations)
        integer(ik), intent(in) :: problem_dimension, population_size, maximum_generations
        real(rk), allocatable :: fitness(:), cov(:,:), chol(:,:), baseline(:,:), baseline_fitness(:), domain_lb(:), domain_ub(:), &
                                 regularization_vector(:), candidates(:,:), candidate_fitness(:), new_cov(:,:)
        real(rk), allocatable, target :: pop1(:,:), pop2(:,:)
        real(rk), pointer :: current_population(:,:), new_population(:,:), dummy_ptr(:,:)
        real(rk) :: mutation_rate, mutation_scale, ft, ft_new, cov_learning_rate, offspring_success, chol_fac, avg_fitness
        integer(ik), allocatable :: selected_pairs_ii(:,:), candidate_sorted_ii(:), improved_average_fitness_ii(:)
        integer(ik) :: generation, elite_ii, i, total_evals, tournament_k, catastrophe_pop_start, catastrophe_count
        logical :: population_ok

        !allocate arrays
        allocate(fitness(population_size), cov(problem_dimension,problem_dimension), chol(problem_dimension,problem_dimension), &
                 domain_lb(problem_dimension), domain_ub(problem_dimension), regularization_vector(problem_dimension), &
                 candidates(problem_dimension,2*population_size), candidate_fitness(2*population_size), &
                 pop1(problem_dimension,population_size), pop2(problem_dimension,population_size), &
                 selected_pairs_ii(2,population_size), candidate_sorted_ii(2*population_size), &
                 new_cov(problem_dimension,problem_dimension), improved_average_fitness_ii(population_size))
!! RASTRIGIN
        domain_lb = -5.12_rk
        domain_ub = 5.12_rk
!!! ROSENBROCK
!        domain_lb = -5.0_rk
!        domain_ub = 10.0_rk

        ! initialize population
        improved_average_fitness_ii = [(i, i=1_ik,population_size)]
        call random_uniform(pop1, size(pop1), minval(domain_lb), maxval(domain_ub))
        current_population => pop1
        new_population => pop2

        ! calculate fitness
        call evaluate_function(current_population, fitness)
        total_evals = size(fitness)
        avg_fitness = sum(fitness)/real(population_size, kind=rk)
        call sort_candidates(fitness, candidate_sorted_ii(1:population_size))
        current_population = current_population(:,candidate_sorted_ii(1:population_size))
        elite_ii = 1
        write(stdout,'(a,f0.6)') 'initial best fitness: ',fitness(elite_ii)
        ft = fitness(elite_ii)


        ! set regularization vector very small, just to avoid numerical collapse
        regularization_vector = (1.0e-10_rk)**2
        ! set mutation rate as 1.0 - 4.0/population_size, enabling high mutation rate for populations 10+
        mutation_rate = 1.0_rk - 4.0_rk/real(population_size, kind=rk)
        ! start mutation scale at 1.0, it will vary depending on generational fitness
        mutation_scale = 1.0_rk
        ! start cov_learning_rate at 0.5, it will vary depending on generational fitness
        cov_learning_rate = 0.05_rk ! 0.05_rk
        ! start tournament_k at 2, it will vary depending on generational fitness
        tournament_k = 2_ik
        ! restart after catastrophe from this position in ranked fitness
        catastrophe_pop_start = population_size/2_ik

        ! perform GA loop
        population_ok = .false.
        catastrophe_count = 0_ik
        do generation=1,maximum_generations

            ! calculate covariance matrix and Cholesky factor for use later
            call covariance(new_cov, current_population(:,improved_average_fitness_ii), reg_vec_opt=regularization_vector)
            if (population_ok) then
                cov = (1.0_rk - cov_learning_rate)*cov + cov_learning_rate*new_cov
            else
                cov = new_cov
            end if
            call cholesky_decomposition(chol, cov)
            
            ! tournament selection
            call perform_selection(current_population, fitness, selected_pairs_ii, k_opt=tournament_k)

            ! store current generation fitness scores
            candidate_fitness(1:population_size) = fitness

            ! fitness-weighted crossover with Gaussian blend
            chol_fac = 0.0_rk
            do i=1_ik,problem_dimension
                chol_fac = chol_fac + cov(i,i)
            end do
            chol_fac = sqrt(real(problem_dimension, kind=rk)/chol_fac)
            call perform_crossover(current_population, selected_pairs_ii, fitness, chol_fac*chol, new_population)
            call apply_constraints(new_population, domain_lb, domain_ub)

            ! Gaussian mutation based on post-crossover population genetic covariance
            call perform_mutation(new_population, mutation_rate, chol, mutation_scale)
            call apply_constraints(new_population, domain_lb, domain_ub)

            ! calculate fitness
            call evaluate_function(new_population, fitness)
            total_evals = total_evals + size(fitness)

            ! store current_population and new_population into candidate_population(2*population_size), then keep top half
            candidates(:,1:population_size) = current_population
            candidates(:,population_size+1:2*population_size) = new_population
            candidate_fitness(population_size+1:2*population_size) = fitness
            call sort_candidates(candidate_fitness, candidate_sorted_ii)
            fitness = candidate_fitness(1:population_size)
            new_population = candidates(:,candidate_sorted_ii(1:population_size))

            ! calculate offspring success = crossover + mutation offspring in population / population size
            offspring_success = real(count(candidate_sorted_ii(1:population_size) > population_size), kind=rk) / &
                                real(population_size, kind=rk)
            if (allocated(improved_average_fitness_ii)) deallocate(improved_average_fitness_ii)
            improved_average_fitness_ii = pack([(i, i=1_ik,population_size)], fitness < avg_fitness)
            avg_fitness = sum(fitness)/real(population_size, kind=rk)
            if (offspring_success >= 0.5_rk) then
                cov_learning_rate = min(1.5_rk*cov_learning_rate, 0.2_rk) ! was 0.1
            else
                cov_learning_rate = max(0.5_rk*cov_learning_rate, 0.01_rk)
            end if

            elite_ii = 1
            ft_new = fitness(elite_ii)
            write(stdout,'(a,i0,3(a,f0.6))') 'generation: ',generation, &
                                             ', best fitness: ',ft_new, &
                                             ', worst fitness: ',fitness(population_size), &
                                             ', offspring success: ',offspring_success

            if (ft_new < ft) then
                mutation_scale = max(0.5_rk*mutation_scale, 0.1_rk) ! was 1.0e-10_rk
            else
                mutation_scale = min(2.0_rk*mutation_scale, 1000.0_rk)
            end if
            ft = ft_new

            if (abs(fitness(1) - fitness(population_size))/fitness(1) > 0.01_rk) then
                population_ok = .true.
            else
                population_ok = .false.
            end if

            if (.not.population_ok) then
                write(*,*) 'entering CATASTROPHE block'
                write(*,*) '  mutation_scale: ',mutation_scale
                write(*,*) '  cov_learning_rate: ',cov_learning_rate
                do concurrent (i=catastrophe_pop_start:population_size)
                    new_population(:,i) = new_population(:,1)
                end do
                mutation_scale = min(mutation_scale*2.0_rk, 1.0_rk)
                call perform_mutation(new_population(:,catastrophe_pop_start:population_size), 1.0_rk, chol, mutation_scale)
                call apply_constraints(new_population(:,catastrophe_pop_start:population_size), domain_lb, domain_ub)
                call evaluate_function(new_population(:,catastrophe_pop_start:population_size), &
                                       fitness(catastrophe_pop_start:population_size))
                total_evals = total_evals + (population_size - catastrophe_pop_start + 1)
                call sort_candidates(fitness, candidate_sorted_ii(1:population_size))
                current_population = new_population(:,candidate_sorted_ii(1:population_size))
                elite_ii = 1
                ft = fitness(elite_ii)
                write(stdout,'(a,i0,2(a,f0.6))') 'CATASTROPHE generation: ',generation, &
                                              ', best fitness: ',ft, &
                                              ', mutation rate: ',mutation_rate
                do i=1,problem_dimension
                    write(*,*) chol(i,i)
                end do
                cov_learning_rate = 0.05_rk ! reset cov_learning_rate to 0.5 for randomized population
                catastrophe_count = catastrophe_count + 1_ik
                if (catastrophe_count > 0_ik) error stop 'death doom loop'
            end if

            ! exit if converged
            if (ft < 1.0e-6) exit

            ! swap population pointers
            dummy_ptr => current_population
            current_population => new_population
            new_population => dummy_ptr
        end do

        ! establish baseline
        allocate(baseline(problem_dimension,total_evals), baseline_fitness(total_evals))
        call random_uniform(baseline, size(baseline), minval(domain_lb), maxval(domain_ub))
        call evaluate_function(baseline, baseline_fitness)
        elite_ii = minloc(baseline_fitness, dim=1)
        write(stdout,'(a,f0.6,a,i0,a)') 'baseline best fitness: ',baseline_fitness(elite_ii),' (',total_evals,' evaluations)'
    end subroutine solve

end module benchmark

program main
use benchmark
implicit none

    integer(ik), parameter :: problem_dimension   = 10
    integer(ik), parameter :: population_size     = problem_dimension*100
    integer(ik), parameter :: maximum_generations = 10

    call solve(problem_dimension, population_size, maximum_generations)

end program main
