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
public :: rastrigin, rosenbrock, griewank, styblinski_tang
public :: real_valued_function
    abstract interface
        pure subroutine real_valued_function(x, fx)
        import rk
        implicit none
            real(rk), intent(in) :: x(:,:)
            real(rk), intent(out) :: fx(:)
        end subroutine real_valued_function
    end interface
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
    pure subroutine griewank(x, fx)
        real(rk), intent(in) :: x(:,:)
        real(rk), intent(out) :: fx(:)
        integer :: i, xdim, nx
        xdim = size(x, dim=1)
        nx = size(x, dim=2)
        do concurrent (i=1:nx)
            fx(i) = 1.0_rk/4000.0_rk*sum(x(:,i)**2, dim=1) - product(cos(x(:,i)/sqrt(real(i, kind=rk)))) + 1.0_rk
        end do
    end subroutine griewank
    pure subroutine styblinski_tang(x, fx)
        real(rk), intent(in) :: x(:,:)
        real(rk), intent(out) :: fx(:)
        integer :: i, xdim, nx
        xdim = size(x, dim=1)
        nx = size(x, dim=2)
        do concurrent (i=1:nx)
            fx(i) = 0.5_rk*sum(x(:,i)**4 - 16.0_rk*x(:,i)**2 + 5.0_rk*x(:,i), dim=1)
        end do
    end subroutine styblinski_tang
end module ga_interface                                                           

module benchmark
use, non_intrinsic :: ga_interface
implicit none
private
public :: stdout, ik, rk, solve, real_valued_function, rastrigin, rosenbrock, griewank, styblinski_tang

    logical, parameter :: print_matrix_enabled = .false., printing = .true.
    integer(ik), parameter :: catastrophe_limit = huge(1_ik)

contains

    impure subroutine print_matrix(mat)
        real(kind=rk), intent(in) :: mat(:,:)
        character(len=32) :: write_format
        integer :: i
        write(write_format, '(a,i0,a)') '(',size(mat,2),'e13.6)'
        do i=1,size(mat,1)
            write(unit=*, fmt=write_format) mat(i,:)
        end do
    end subroutine print_matrix

    impure subroutine solve(evaluate_function, target_value, domain_lb, domain_ub)
        procedure(real_valued_function) :: evaluate_function
        real(rk), intent(in) :: target_value, domain_lb(:), domain_ub(:)
        real(rk), allocatable :: fitness(:), cov(:,:), chol(:,:),  &
                                 regularization_vector(:), candidates(:,:), candidate_fitness(:), new_cov(:,:), cov_weights(:)
        real(rk), allocatable, target :: pop1(:,:), pop2(:,:)
        real(rk), pointer :: current_population(:,:), new_population(:,:), dummy_ptr(:,:)
        real(rk) :: mutation_rate, mutation_scale, cov_learning_rate, offspring_success, chol_fac, &
                    current_fitness_stats(4), new_fitness_stats(4), & ! best, median, worst, average
                    mutation_scale0, mutation_scale_min, mutation_scale_max, &
                    cov_learning_rate0, cov_learning_rate_min, cov_learning_rate_max
        integer(ik), allocatable :: selected_pairs_ii(:,:), candidate_sorted_ii(:)
        integer(ik) :: generation, i, total_evals, tournament_k, catastrophe_pop_start, catastrophe_count, &
                       problem_dimension, population_size, maximum_generations
        logical :: population_ok

        problem_dimension = size(domain_lb)
        population_size = 10_ik + 10_ik*problem_dimension
        maximum_generations = int(1000000.0_rk/real(population_size, kind=rk), kind=ik)

        !allocate arrays
        allocate(fitness(population_size), cov(problem_dimension,problem_dimension), chol(problem_dimension,problem_dimension), &
                 regularization_vector(problem_dimension), &
                 candidates(problem_dimension,2*population_size), candidate_fitness(2*population_size), &
                 pop1(problem_dimension,population_size), pop2(problem_dimension,population_size), &
                 selected_pairs_ii(2,population_size), candidate_sorted_ii(2*population_size), &
                 new_cov(problem_dimension,problem_dimension), cov_weights(population_size))

        ! initialize population
        call random_uniform(pop1, size(pop1), minval(domain_lb), maxval(domain_ub))
        current_population => pop1
        new_population => pop2

        ! calculate fitness
        call evaluate_function(current_population, fitness)
        fitness = abs(target_value - fitness)
        total_evals = size(fitness)
        call sort_candidates(fitness, candidate_sorted_ii(1:population_size))
        current_population = current_population(:,candidate_sorted_ii(1:population_size))
        current_fitness_stats = [fitness(1), &
                                 fitness(population_size/2_ik), &
                                 fitness(population_size), &
                                 sum(fitness)/real(population_size, kind=rk)]
        if (printing) then
            write(stdout,'(a,f0.6)') 'initial best fitness: ',current_fitness_stats(1)
        end if

        ! set regularization vector very small, just to avoid numerical collapse
        regularization_vector = (1.0e-10_rk)**2

        ! set mutation rate as 1.0 - 4.0/population_size, enabling high mutation rate for populations 10+
        mutation_rate = 0.1_rk ! 1.0_rk - 4.0_rk/real(population_size, kind=rk)

        ! start mutation scale at 1.0, it will vary depending on generational fitness
        mutation_scale0 = (1.0_rk/6.0_rk)/sqrt(real(problem_dimension, kind=rk))
        mutation_scale = mutation_scale0
        mutation_scale_min = 0.5_rk*mutation_scale0
        mutation_scale_max = 3.0_rk*mutation_scale0
        
        ! start cov_learning_rate at 0.5, it will vary depending on generational fitness
        cov_learning_rate0 = 0.05_rk
        cov_learning_rate = cov_learning_rate0
        cov_learning_rate_min = 0.01_rk
        cov_learning_rate_max = 0.5_rk

        ! start tournament_k at 2, it will vary depending on generational fitness
        tournament_k = 2_ik
        ! restart after catastrophe from this position in ranked fitness
        catastrophe_pop_start = population_size/2_ik

        ! perform GA loop
        population_ok = .false.
        catastrophe_count = 0_ik
        do generation=1,maximum_generations

            ! calculate covariance matrix and Cholesky factor for use later
            cov_weights = log(real(population_size, kind=rk) + 1.0_rk) - log([(real(i, kind=rk), i=1_ik,population_size)])
            cov_weights = cov_weights/sum(cov_weights)
            where (cov_weights < 0.001_rk) cov_weights = 0.0_rk
            cov_weights = cov_weights/sum(cov_weights)
            call covariance(new_cov, current_population, reg_vec_opt=regularization_vector, weights_opt=cov_weights)
            if (print_matrix_enabled) call print_matrix(cov)
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
            chol = 0.0_rk
            do concurrent (i=1_ik:problem_dimension)
                chol(i,i) = domain_ub(i) - domain_lb(i)
            end do
            call perform_mutation(new_population, mutation_rate, chol, mutation_scale)
            call apply_constraints(new_population, domain_lb, domain_ub)

            ! calculate fitness
            call evaluate_function(new_population, fitness)
            fitness = abs(target_value - fitness)
            total_evals = total_evals + size(fitness)

            ! store current_population and new_population into candidate_population(2*population_size), then keep top half
            candidates(:,1:population_size) = current_population
            candidates(:,population_size+1:2*population_size) = new_population
            candidate_fitness(population_size+1:2*population_size) = fitness
            call sort_candidates(candidate_fitness, candidate_sorted_ii)
            fitness = candidate_fitness(1:population_size)
            new_population = candidates(:,candidate_sorted_ii(1:population_size))

            new_fitness_stats = [fitness(1), &
                                 fitness(population_size/2_ik), &
                                 fitness(population_size), &
                                 sum(fitness)/real(population_size, kind=rk)]

            ! calculate offspring success = crossover + mutation offspring in population / population size
            offspring_success = real(count(candidate_sorted_ii(1:population_size) > population_size), kind=rk) / &
                                real(population_size, kind=rk)

            ! ~18.4% is optimal for some valley function, ~20% standard to balance convergence speed while not overshooting
            if (offspring_success >= 0.2_rk) then
                cov_learning_rate = min(1.5_rk*cov_learning_rate, cov_learning_rate_max)
                tournament_k = max(tournament_k - 1_ik, 2_ik)
            else
                cov_learning_rate = max(0.5_rk*cov_learning_rate, cov_learning_rate_min)
            end if

            ! spread from median fitness to best fitness captures top half of distribution, if it is collapsing, we need more noise
            if ((new_fitness_stats(1) < current_fitness_stats(1)) .or. &
                ((new_fitness_stats(2) - new_fitness_stats(1)) > (current_fitness_stats(2) - current_fitness_stats(1)))) then
                mutation_scale = max(0.5_rk*mutation_scale, mutation_scale_min)
            else
                mutation_scale = min(1.5_rk*mutation_scale, mutation_scale_max)
            end if

            if (printing) then
                write(stdout,'(a,i0,4(a,f0.6))') '  generation: ',generation, &
                                                 ', best fitness: ',new_fitness_stats(1), &
                                                 ', offspring success: ',offspring_success, &
                                                 ', cov learning rate: ',cov_learning_rate, &
                                                 ', mutation scale: ',mutation_scale
            end if

            if ((mutation_scale >= 0.7) .and. (new_fitness_stats(1) >= current_fitness_stats(1))) then
                population_ok = .false.
            else
                population_ok = .true.
                catastrophe_count = 0_ik
            end if

            if (.not.population_ok) then
                if (printing) then
                    write(stdout,'(a)') '    CATASTROPHE - doubling covariance'
                    write(stdout,'(a,4e13.6)') '      previous population fitness min/med/max/avg: ',current_fitness_stats
                    write(stdout,'(a,4e13.6)') '           new population fitness min/med/max/avg: ',new_fitness_stats
                end if

                cov = 2.0_rk*cov ! double covariance to help crossover search
                mutation_scale = mutation_scale0 ! reset mutation_scale to 0.5
                cov_learning_rate = cov_learning_rate_min ! reset cov_learning_rate to 0.01, prevent immediate collapse
                tournament_k = tournament_k + 1_ik ! increase tournament size to increase selection pressure
                do concurrent (i=2_ik:population_size)
                    new_population(:,i) = new_population(:,1_ik)
                end do
                call perform_mutation(new_population(:,2:population_size), 1.0_rk, chol, mutation_scale) ! copy + mutate elite 1

                catastrophe_count = catastrophe_count + 1_ik
                if (catastrophe_count > catastrophe_limit) error stop 'death doom loop'
            end if

            ! exit if converged
            current_fitness_stats = new_fitness_stats
            if (current_fitness_stats(1) < 1.0e-6) then
                exit
            end if

            ! swap population pointers
            dummy_ptr => current_population
            current_population => new_population
            new_population => dummy_ptr
        end do
        write(stdout,*) '  value: ',current_fitness_stats(1),', evaluations: ',total_evals
    end subroutine solve

end module benchmark

program main
use benchmark
implicit none

    integer(ik), parameter :: k_test_functions = 4
    integer(ik), parameter :: d_dimension_list(*) = [2_ik, 10_ik, 20_ik, 100_ik, 200_ik]

    procedure(real_valued_function), pointer :: test_function
    character(len=:), allocatable :: fname
    real(rk) :: target_value
    real(rk), allocatable :: domain_lb(:), domain_ub(:)
    integer(ik) :: k, d_ii, d

    do d_ii=1_ik,size(d_dimension_list)
        d = d_dimension_list(d_ii)
        if (allocated(domain_lb)) deallocate(domain_lb)
        if (allocated(domain_ub)) deallocate(domain_ub)
        allocate(domain_lb(d), domain_ub(d))
        do k=1_ik,k_test_functions
            select case (k)
                case (1)
                    test_function => rastrigin
                    fname = 'Rastrigin'
                    target_value = 0.0_rk
                    domain_lb = -5.12_rk
                    domain_ub = 5.12_rk
                case (2)
                    test_function => rosenbrock
                    fname = 'Rosenbrock'
                    target_value = 0.0_rk
                    domain_lb = -5.0_rk
                    domain_ub = 10.0_rk
                case (3)
                    test_function => griewank
                    fname = 'Griewank'
                    target_value = 0.0_rk
                    domain_lb = 100.0_rk
                    domain_ub = 100.0_rk
                case (4)
                    test_function => styblinski_tang
                    fname = 'Styblinski-Tang'
                    target_value = -39.166165_rk*real(d, kind=rk)
                    domain_lb = -5.0_rk
                    domain_ub = 5.0_rk
                case default
                    write(stdout,'(a)') 'valid test functions are:'
                    write(stdout,'(a)') ' 1. Rastrigin'
                    write(stdout,'(a)') ' 2. Rosenbrock'
                    write(stdout,'(a)') ' 3. Griewank'
                    write(stdout,'(a)') ' 4. Styblinski-Tang'
                    error stop 'only 4 test functions have been implemented'
            end select
            write(stdout,'(a,f0.6,a,i0,a,f0.2,a,f0.2,a)') &
                 fname//' looking for ',target_value,' on ',d,' dimensions [',minval(domain_lb),', ',maxval(domain_ub),']'
            call solve(test_function, target_value, domain_lb, domain_ub)
        end do
    end do

end program main
