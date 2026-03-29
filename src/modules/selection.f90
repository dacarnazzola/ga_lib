module selection
use, non_intrinsic :: kinds, only: i32, i64, sp
use, non_intrinsic :: system, only: debug_error_condition
use, non_intrinsic :: random, only: random_uniform_i32
implicit none
private
public :: selection_uniform_random, selection_tournament

    integer(kind=i32), parameter :: DEFAULT_TOURNAMENT_SIZE = 2

contains

    impure subroutine selection_uniform_random(population, selected_pairs_ii)
        real(kind=sp), intent(in) :: population(:,:)
        integer(kind=i32), intent(out) :: selected_pairs_ii(2,size(population, dim=2, kind=i64))
        call debug_error_condition(size(population,dim=2,kind=i64) > huge(1_i32), &
                                   'SELECTION::SELECTION_UNIFORM_RANDOM population is too large for i32 storage')
        call random_uniform_i32(selected_pairs_ii, size(selected_pairs_ii, kind=i32), 1_i32, size(population, dim=2, kind=i32))
    end subroutine selection_uniform_random

    impure subroutine selection_tournament(population, fitness, selected_pairs_ii, k_opt)
        real(kind=sp), intent(in) :: population(:,:), fitness(size(population, dim=2, kind=i64))
        integer(kind=i32), intent(out) :: selected_pairs_ii(2,size(population, dim=2, kind=i64))
        integer(kind=i32), intent(in), optional :: k_opt
        integer(kind=i32) :: k, i, t1_start, t1_end, t2_start, t2_end
        integer(kind=i32), allocatable :: tournaments(:,:)
        real(kind=sp), allocatable :: tournament_fitness(:)
        call debug_error_condition(size(population,dim=2,kind=i64) > huge(1_i32), &
                                   'SELECTION::SELECTION_TOURNAMENT population is too large for i32 storage')
        call debug_error_condition(size(population,dim=2) /= size(fitness), &
                                   'SELECTION::SELECTION_TOURNAMENT population size does not match fitness size')
        if (present(k_opt)) then
            call debug_error_condition(k_opt < 1_i32, &
                                       'SELECTON::SELECTION_TOURNAMENT k must be >= 1')
            k = k_opt
        else
            k = DEFAULT_TOURNAMENT_SIZE
        end if
        t1_start = 1_i32
        t1_end = k
        t2_start = k + 1_i32
        t2_end = 2_i32*k
        allocate(tournaments(2_i32*k,size(population,dim=2,kind=i32)))
        call random_uniform_i32(tournaments, size(tournaments, kind=i32), 1_i32, size(population, dim=2, kind=i32))
        allocate(tournament_fitness(2_i32*k))
        do i=1_i32,size(population, dim=2, kind=i32)
            tournament_fitness = fitness(tournaments(:,i))
            selected_pairs_ii(1,i) = tournaments(minloc(tournament_fitness(t1_start:t1_end), dim=1),i)
            selected_pairs_ii(2,i) = tournaments(minloc(tournament_fitness(t2_start:t2_end), dim=1)+k,i)
        end do
    end subroutine selection_tournament

end module selection
