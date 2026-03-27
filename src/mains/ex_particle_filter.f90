module pf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: pi_dp, twopi_dp, deg2rad_dp, rad2deg_dp, eps_dp
use, non_intrinsic :: random, only: random_normal, random_uniform, random_log_uniform
use, non_intrinsic :: statistics, only: dsum, avg, std, cumsum, cov, normalize, geomean
use, non_intrinsic :: vector_math, only: vdot
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: sorting, only: sort
implicit none
private

    logical, parameter :: debug      = .false.
    logical, parameter :: debug_run  = .false.
    logical, parameter :: debug_init = .false.
    logical :: ballistic_particles   = .false.
    integer :: vt_initialization_mode = -1

    integer, parameter :: max_iters = 100
    
    real(dp), parameter :: global_minimum_range = 1.0_dp
    real(dp), parameter :: nmi2ft = 1852.0_dp*100.0_dp/2.54_dp/12.0_dp
    real(dp), parameter :: ft2nmi = 12.0_dp*2.54_dp/100.0_dp/1852.0_dp
    real(dp), parameter :: g = 32.2_dp

    interface
        impure subroutine resample_cartesian_particles(n, cartesian_particles, weights, jerk_sig)
        import dp
        implicit none
            integer, intent(in) :: n
            real(dp), intent(inout) :: cartesian_particles(:,:), weights(:), jerk_sig(:)
        end
        pure subroutine calculate_weights(tgt_pol, meas_sig, n, polar_particles, weights)
        import dp
        implicit none
            real(dp), intent(in) :: tgt_pol(:), meas_sig(:)
            integer, intent(in) :: n
            real(dp), intent(in) :: polar_particles(:,:)
            real(dp), intent(inout) :: weights(:)
        end
    end interface

    public :: debug, debug_run, debug_init, dp, nmi2ft, ft2nmi, deg2rad_dp, rad2deg_dp, g, pi_dp, twopi_dp, &
              perfect_cart2pol, generate_measurements, initialize_particles, convert_particles_cart2pol, &
              generate_state_estimate, calculate_weights_gaussian, calculate_weights_student_t, regularize_covariance, &
              calculate_kernel_bandwidth, apply_kernel_smoothing, calculate_weights, &
              resample_cartesian_particles, resample_systematic, resample_multinomial, resample_residual, resample_stratified, &
              fly_constant_acceleration, propagate_particles, apply_constraints, resample_jerk, dump_particles, &
              rmse, neff, &
              avg, std, nearly, geomean, sort, &
              ballistic_particles, vt_initialization_mode

contains

    pure subroutine perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
        real(dp), intent(in) :: obs_cart(:), tgt_cart(size(obs_cart))
        real(dp), intent(out) :: tgt_pol(3)
        call cart2pol_inner(obs_cart, tgt_cart, 0.0_dp, 0.0_dp, 0.0_dp, tgt_pol)
    end

    pure subroutine cart2pol_inner(obs_cart, tgt_cart, err_rng, err_ang, err_rngrt, tgt_pol)
        real(dp), intent(in) :: obs_cart(:), tgt_cart(size(obs_cart)), err_rng, err_ang, err_rngrt
        real(dp), intent(out) :: tgt_pol(3)
        real(dp) :: rng, loshat(2), losv(2)
        rng = sqrt((tgt_cart(1) - obs_cart(1))**2 + (tgt_cart(2) - obs_cart(2))**2)
        tgt_pol(1) = max(global_minimum_range, rng + err_rng)
        loshat = (tgt_cart(1:2) - obs_cart(1:2))/rng
        tgt_pol(2) = mod(atan2(loshat(2), loshat(1)) + err_ang + pi_dp, twopi_dp) - pi_dp
        losv = tgt_cart(3:4) - obs_cart(3:4)
        tgt_pol(3) = losv(1)*loshat(1) + losv(2)*loshat(2) + err_rngrt
    end

    impure subroutine generate_measurements(obs_cart, tgt_cart, meas_sig, no_meas_max_rng, no_meas_max_spd, n, polar_measurements)
        real(dp), intent(in) :: obs_cart(:), tgt_cart(size(obs_cart)), meas_sig(3), no_meas_max_rng, no_meas_max_spd
        integer, intent(in) :: n
        real(dp), intent(out) :: polar_measurements(3,n)
        real(dp) :: err_rng(n), err_ang(n), err_rngrt(n), true_rng, true_ang, losv(2), loshat(2), true_rngrt
        integer :: i
        if (meas_sig(1) > 0) then
            call random_normal(err_rng, 0.0_dp, meas_sig(1))
        else
            true_rng = max(1.0e-6_dp, sqrt((tgt_cart(1) - obs_cart(1))**2 + (tgt_cart(2) - obs_cart(2))**2))
            call random_uniform(err_rng, global_minimum_range - true_rng, no_meas_max_rng - true_rng)
        end if
        if (meas_sig(2) > 0) then
            call random_normal(err_ang, 0.0_dp, meas_sig(2))
        else
            true_ang = atan2(tgt_cart(2) - obs_cart(2), tgt_cart(1) - obs_cart(1))
            call random_uniform(err_ang, -pi_dp - true_ang, pi_dp - true_ang)
        end if
        if (meas_sig(3) > 0) then
            call random_normal(err_rngrt, 0.0_dp, meas_sig(3))
        else
            losv = tgt_cart(3:4) - obs_cart(3:4)
            loshat = (tgt_cart(1:2) - obs_cart(1:2)) / &
                     max(1.0e-6_dp, sqrt((tgt_cart(1) - obs_cart(1))**2 + (tgt_cart(2) - obs_cart(2))**2))
            true_rngrt = losv(1)*loshat(1) + losv(2)*loshat(2)
            call random_uniform(err_rngrt, -no_meas_max_spd - true_rngrt, no_meas_max_spd - true_rngrt)
        end if
        do concurrent (i=1:n)
            call cart2pol_inner(obs_cart, tgt_cart, err_rng(i), err_ang(i), err_rngrt(i), polar_measurements(:,i))
        end do
    end

    impure subroutine resample_measurements(meas, meas_sig, no_meas_max_rng, no_meas_max_spd, n, polar_measurements)
        real(dp), intent(in) :: meas(3), meas_sig(3), no_meas_max_rng, no_meas_max_spd
        integer, intent(in) :: n
        real(dp), intent(out) :: polar_measurements(3,n)
        integer :: i, iter
        if (meas_sig(1) > 0) then
            call random_normal(polar_measurements(1,:), meas(1), meas_sig(1))
        else
            call random_uniform(polar_measurements(1,:), global_minimum_range, no_meas_max_rng)
        end if
        if (meas_sig(2) > 0) then
            call random_normal(polar_measurements(2,:), meas(2), meas_sig(2))
        else
            call random_uniform(polar_measurements(2,:), -pi_dp, pi_dp)
        end if
        if (ballistic_particles) then
            do i=1,n
                iter = 0
                loop_reroll: do while (polar_measurements(2,i) < 0.0_dp)
                    call random_normal(polar_measurements(2:2,i), meas(2), meas_sig(2))
                    iter = iter + 1
                    if (iter > max_iters) exit loop_reroll
                end do loop_reroll
            end do
        end if
        if (meas_sig(3) > 0) then
            call random_normal(polar_measurements(3,:), meas(3), meas_sig(3))
        else
            call random_uniform(polar_measurements(3,:), -no_meas_max_spd, no_meas_max_spd)
        end if
    end

    impure subroutine resample_jerk(jerk_sig, lo, hi)
        real(dp), intent(out) :: jerk_sig(:)
        real(dp), intent(in), optional :: lo, hi
        real(dp) :: jerk_sig_lo, jerk_sig_hi
        if (present(lo)) then
            jerk_sig_lo = lo
        else
            jerk_sig_lo = 10.0_dp**(-3.25_dp) !! yields p10 value 0.001, p50 value 0.01, p90 value 0.1 with both defaults
        end if
        if (present(hi)) then
            jerk_sig_hi = hi
        else
            jerk_sig_hi = 10.0_dp**(-0.75_dp) !! yields p10 value 0.001, p50 value 0.01, p90 value 0.1 with both defaults
        end if
        call random_log_uniform(jerk_sig, jerk_sig_lo, jerk_sig_hi)
        jerk_sig = jerk_sig*g
    end

    impure subroutine initialize_particles(obs_cart, meas, meas_sig, tgt_max_rng, tgt_max_spd, n, &
                                           cartesian_particles, weights, jerk_sig)
        real(dp), intent(in) :: obs_cart(:), meas(3), meas_sig(3), tgt_max_rng, tgt_max_spd
        integer, intent(in) :: n
        real(dp), intent(out) :: cartesian_particles(size(obs_cart),n), weights(n), jerk_sig(n)
        real(dp) :: polar_measurements(3,n), tgt_spd_scale(n), ax(n), ay(n)
        integer :: i, iter
        call resample_measurements(meas, meas_sig, tgt_max_rng, tgt_max_spd, n, polar_measurements)
        select case (vt_initialization_mode)
            case (1)
                call random_uniform(tgt_spd_scale, 0.0_dp, 1.0_dp)
            case (2)
                call random_log_uniform(tgt_spd_scale, 0.25_dp, 1.0_dp) !! p50 value of 0.5
            case (3)
                call random_normal(tgt_spd_scale, 0.0_dp, 1.0_dp)
                do i=1,n
                    iter = 0
                    loop_reroll: do while ((tgt_spd_scale(i) > 1.0_dp) .or. (tgt_spd_scale(i) < -1.0_dp))
                        call random_normal(tgt_spd_scale(i:i), 0.0_dp, 1.0_dp)
                        iter = iter + 1
                        if (iter > max_iters) exit loop_reroll
                    end do loop_reroll
                end do
                tgt_spd_scale = abs(tgt_spd_scale)
            case (4)
                call random_normal(tgt_spd_scale, 0.0_dp, 1.0_dp/2.0_dp)
                do i=1,n
                    iter = 0
                    loop_reroll2: do while ((tgt_spd_scale(i) > 1.0_dp) .or. (tgt_spd_scale(i) < -1.0_dp))
                        call random_normal(tgt_spd_scale(i:i), 0.0_dp, 1.0_dp/2.0_dp) 
                        iter = iter + 1
                        if (iter > max_iters) exit loop_reroll2
                    end do loop_reroll2
                end do
                tgt_spd_scale = abs(tgt_spd_scale)
            case (5)
                call random_normal(tgt_spd_scale, 0.0_dp, 1.0_dp/3.0_dp)
                do i=1,n
                    loop_reroll3: do while ((tgt_spd_scale(i) > 1.0_dp) .or. (tgt_spd_scale(i) < -1.0_dp))
                        call random_normal(tgt_spd_scale(i:i), 0.0_dp, 1.0_dp/3.0_dp)
                        iter = iter + 1
                        if (iter > max_iters) exit loop_reroll3
                    end do loop_reroll3
                end do
                tgt_spd_scale = abs(tgt_spd_scale)
            case default
                error stop 'unimplemented Vtangential initialization mode'
        end select
        call resample_jerk(jerk_sig)
        if (ballistic_particles) then
            ax = 0.0_dp
            ay = -1.0_dp*g
        else
            call random_normal(ax, 0.0_dp, 1.0_dp*g)
            call random_normal(ay, 0.0_dp, 1.0_dp*g)
        end if
        do concurrent (i=1:n)
            call pol2cart_inner(obs_cart, tgt_max_rng, tgt_max_spd, (-1)**mod(i,2)*tgt_spd_scale(i), &
                                polar_measurements(:,i), cartesian_particles(:,i))
            if (size(obs_cart) == 6) then
                cartesian_particles(5,i) = ax(i)
                cartesian_particles(6,i) = ay(i)
            end if
            if (ballistic_particles) then
                cartesian_particles(2,i) = max(0.0_dp, cartesian_particles(2,i))
                cartesian_particles(4,i) = abs(cartesian_particles(4,i))
            end if
        end do
        weights = 1.0_dp/real(n,dp)
    end

    pure subroutine pol2cart_inner(obs_cart, tgt_max_rng, tgt_max_spd, tgt_spd_scale, tgt_pol, tgt_cart)
        real(dp), intent(in) :: obs_cart(:), tgt_max_rng, tgt_max_spd, tgt_spd_scale, tgt_pol(3)
        real(dp), intent(out) :: tgt_cart(size(obs_cart))
        real(dp) :: cos_ang, sin_ang, tgt_rng, tgt_min_spd, tgt_spd, tgt_vt_mag, tgt_final_spd_scale
        integer :: tgt_vt_dir
        cos_ang = cos(tgt_pol(2))
        sin_ang = sin(tgt_pol(2))
        tgt_rng = min(tgt_max_rng, tgt_pol(1))
        tgt_cart(1) = obs_cart(1) + tgt_rng*cos_ang
        tgt_cart(2) = obs_cart(2) + tgt_rng*sin_ang
        tgt_cart(3) = obs_cart(3) + tgt_pol(3)*cos_ang
        tgt_cart(4) = obs_cart(4) + tgt_pol(3)*sin_ang
        tgt_min_spd = min(tgt_max_spd, sqrt(tgt_cart(3)**2 + tgt_cart(4)**2))
        tgt_spd = tgt_spd_scale*(tgt_max_spd - tgt_min_spd) + tgt_min_spd
        tgt_vt_mag = sqrt(max(0.0_dp, tgt_spd**2 - tgt_min_spd**2))
        tgt_vt_dir = (-1)**mod(floor(tgt_spd_scale*100000000.0_dp), 2)
        tgt_cart(3) = tgt_cart(3) - tgt_vt_dir*tgt_vt_mag*sin_ang
        tgt_cart(4) = tgt_cart(4) + tgt_vt_dir*tgt_vt_mag*cos_ang
        tgt_final_spd_scale = max(1.0_dp, sqrt(tgt_cart(3)**2 + tgt_cart(4)**2)/tgt_spd)
        tgt_cart(3:4) = tgt_cart(3:4)/tgt_final_spd_scale
    end

    pure subroutine convert_particles_cart2pol(obs_cart, n, cartesian_particles, polar_particles)
        real(dp), intent(in) :: obs_cart(:) 
        integer, intent(in) :: n
        real(dp), intent(in) :: cartesian_particles(size(obs_cart),n)
        real(dp), intent(out) :: polar_particles(3,n)
        integer :: i
        do concurrent (i=1:n)
            call perfect_cart2pol(obs_cart, cartesian_particles(:,i), polar_particles(:,i))
        end do
    end

    pure subroutine generate_state_estimate(state_estimate, n, particles, weights)
        real(dp), intent(out) :: state_estimate(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: particles(size(state_estimate),n), weights(n)
        integer :: i
        do concurrent (i=1:size(state_estimate))
            state_estimate(i) = vdot(particles(i,:), weights)
        end do
    end

    pure subroutine calculate_weights_gaussian(tgt_pol, meas_sig, n, polar_particles, weights)
        real(dp), intent(in) :: tgt_pol(:), meas_sig(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: polar_particles(:,:)
        real(dp), intent(inout) :: weights(:)
        real(dp) :: err_fac(3), log_weights(n)
        integer :: i
        call debug_error_condition(size(tgt_pol) /= 3, 'assumed polar: range, angle, range-rate')
        call debug_error_condition(size(tgt_pol) /= size(meas_sig), 'tgt_pol and meas_sig size must match')
        call debug_error_condition(size(tgt_pol) /= size(polar_particles,dim=1), 'tgt_pol and polar_particle ndims must match')
        call debug_error_condition(size(weights) /= size(polar_particles,dim=2), 'size of weights and polar_particles nsamples')
        do concurrent (i=1:n)
            err_fac = 0.0_dp
            if (meas_sig(1) > 0) err_fac(1) = (polar_particles(1,i) - tgt_pol(1))**2/meas_sig(1)**2
            if (meas_sig(2) > 0) err_fac(2) = (mod(polar_particles(2,i) - tgt_pol(2) + pi_dp, twopi_dp) - pi_dp)**2/meas_sig(2)**2
            if (meas_sig(3) > 0) err_fac(3) = (polar_particles(3,i) - tgt_pol(3))**2/meas_sig(3)**2
            log_weights(i) = -0.5_dp*(err_fac(1) + err_fac(2) + err_fac(3)) !! log domain
        end do
        log_weights = log_weights - maxval(log_weights) !! shift log domain weights by maxmimum log-weight
        weights = weights*exp(log_weights) ! + eps_dp !! convert back to non-log domain, consider adding eps_dp if this is still collapsing
        weights = weights/dsum(weights)
    end

    pure subroutine calculate_weights_student_t(tgt_pol, meas_sig, n, polar_particles, weights)
        real(dp), intent(in) :: tgt_pol(:), meas_sig(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: polar_particles(:,:)
        real(dp), intent(inout) :: weights(:)
        real(dp) :: err_fac(3), log_weights(n), v(n)
        integer :: i
        call debug_error_condition(size(tgt_pol) /= 3, 'assumed polar: range, angle, range-rate')
        call debug_error_condition(size(tgt_pol) /= size(meas_sig), 'tgt_pol and meas_sig size must match')
        call debug_error_condition(size(tgt_pol) /= size(polar_particles,dim=1), 'tgt_pol and polar_particle ndims must match')
        call debug_error_condition(size(weights) /= size(polar_particles,dim=2), 'size of weights and polar_particles nsamples')
        call debug_error_condition(minval(meas_sig) <= 0.0_dp, 'student_t weights not set up for missing dimensions')
        do concurrent (i=1:n)
            v(i) = 1.0_dp + (meas_sig(1)**2)/((polar_particles(1,i)**2)*(meas_sig(2)**2))
            err_fac(1) = log(1.0_dp + (polar_particles(1,i) - tgt_pol(1))**2/(v(i)*meas_sig(1)**2))
            err_fac(2) = log(1.0_dp + (mod(polar_particles(2,i) - tgt_pol(2) + pi_dp, twopi_dp) - pi_dp)**2/(v(i)*meas_sig(2)**2))
            err_fac(3) = log(1.0_dp + (polar_particles(3,i) - tgt_pol(3))**2/(v(i)*meas_sig(3)**2))
            log_weights(i) = (-(v(i) + 1.0_dp)/2.0_dp)*(err_fac(1) + err_fac(2) + err_fac(3)) !! log domain
        end do
        log_weights = log_weights - maxval(log_weights) !! shift log domain weights by maxmimum log-weight
        weights = weights*exp(log_weights) ! + eps_dp !! convert back to non-log domain, consider adding eps_dp if this is still collapsing
        weights = weights/dsum(weights)
!        if (debug) then
!            write(*,*) 'min v: ',minval(v),', max v: ',maxval(v)
!            write(*,*) 'min log_weight: ',minval(log_weights),', max log_weight: ',maxval(log_weights)
!            write(*,*) 'min weight: ',minval(weights),', max weight: ',maxval(weights)
!        end if
    end

    impure subroutine resample_systematic(n, cartesian_particles, weights, jerk_sig)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:), jerk_sig(:)
        real(dp) :: u(1), inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2)), &
                    new_jerk_sig(size(jerk_sig))
        integer :: i, j
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        call debug_error_condition(size(jerk_sig) /= n, 'size of jerk_sig array needs to match n')
        call cumsum(weights)
        weights(n) = 1.0_dp
        inv_n = 1.0_dp/real(n, kind=dp)
        call random_uniform(u, 0.0_dp, inv_n)
        j = 1
        do i=1,n
            do while (u(1) > weights(j))
                j = j + 1
            end do
            new_particles(:,i) = cartesian_particles(:,j)
            new_jerk_sig(i) = jerk_sig(j)
            u = u + inv_n
        end do
        cartesian_particles = new_particles
        jerk_sig = new_jerk_sig
        weights = inv_n
    end

    impure subroutine resample_multinomial(n, cartesian_particles, weights, jerk_sig)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:), jerk_sig(:)
        real(dp) :: u(n), inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2)), &
                    new_jerk_sig(size(jerk_sig))
        integer :: i, j
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        call debug_error_condition(size(jerk_sig) /= n, 'size of jerk_sig array needs to match n')
        call cumsum(weights)
        weights(n) = 1.0_dp
        call random_uniform(u, 0.0_dp, 1.0_dp)
        do i=1,n
            j = 1
            do while (u(i) > weights(j))
                j = j + 1
            end do
            new_particles(:,i) = cartesian_particles(:,j)
            new_jerk_sig(i) = jerk_sig(j)
        end do
        cartesian_particles = new_particles
        jerk_sig = new_jerk_sig
        inv_n = 1.0_dp/real(n, kind=dp)
        weights = inv_n
    end

    impure subroutine resample_residual(n, cartesian_particles, weights, jerk_sig)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:), jerk_sig(:)
        real(dp) :: inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2)), resamples(n), u(n), &
                    new_jerk_sig(size(jerk_sig))
        integer :: i, resamples_int, j, k
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        call debug_error_condition(size(jerk_sig) /= n, 'size of jerk_sig array needs to match n')
        resamples = weights*n
        j = 1
        do i=1,n
            resamples_int = floor(resamples(i))
            do k=1,resamples_int
                new_particles(:,j) = cartesian_particles(:,i)
                new_jerk_sig(j) = jerk_sig(i)
                j = j + 1
            end do
        end do
        resamples = resamples - real(floor(resamples), kind=dp)
        if (dsum(resamples) > 0.0_dp) then
            call normalize(resamples)
            call cumsum(resamples)
            resamples(n) = 1.0_dp
            call random_uniform(u(1:n-j+1), 0.0_dp, 1.0_dp)
            do i=1,n-j+1
                k = 1
                do while (u(i) > resamples(k))
                    k = k + 1
                end do
                new_particles(:,i+j-1) = cartesian_particles(:,k)
                new_jerk_sig(i+j-1) = jerk_sig(k)
            end do
        end if
        cartesian_particles = new_particles
        jerk_sig = new_jerk_sig
        inv_n = 1.0_dp/real(n, kind=dp)
        weights = inv_n
    end

    impure subroutine resample_stratified(n, cartesian_particles, weights, jerk_sig)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:), jerk_sig(:)
        real(dp) :: u(n), inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2)), &
                    new_jerk_sig(size(jerk_sig))
        integer :: i, j
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        inv_n = 1.0_dp/real(n, kind=dp)
        call random_uniform(u, 0.0_dp, 1.0_dp)
        do i=1,n
            u(i) = (u(i) + real(i-1,dp))*inv_n
        end do
        call cumsum(weights)
        weights(n) = 1.0_dp
        j = 1
        do i=1,n
            do while (u(i) > weights(j))
                j = j + 1
            end do
            new_particles(:,i) = cartesian_particles(:,j)
            new_jerk_sig(i) = jerk_sig(j)
        end do
        cartesian_particles = new_particles
        jerk_sig = new_jerk_sig
        weights = inv_n
    end

    impure subroutine regularize_covariance(z_scale, n, particles, weights, jerk_sig)
        real(dp), intent(in) :: z_scale
        integer, intent(in) :: n
        real(dp), intent(inout) :: particles(:,:), weights(:), jerk_sig(:)
        real(dp) :: particles_cov(size(particles,dim=1),size(particles,dim=1)), &
                    mu_zero(size(particles,dim=1)), z(size(particles,dim=1),size(particles,dim=2)), &
                    new_particles(size(particles,dim=1),size(particles,dim=2))
        integer :: i, particles_ii(n)
        call debug_error_condition(size(particles, dim=2) /= n, 'must have n particles')
        call debug_error_condition(size(weights) /= n, 'size of weights must match n')
        call debug_error_condition(size(jerk_sig) /= n, 'size of jerk_sig must match n')
        call cov(particles, weights, particles_cov)
        mu_zero = 0.0_dp
        call random_normal(z, mu_zero, particles_cov)
        call random_uniform(particles_ii, 1, n)
        do concurrent (i=1:n)
            new_particles(:,i) = particles(:,particles_ii(i)) + z_scale*z(:,i)
        end do
        particles = new_particles
        weights = 1.0_dp/real(n, kind=dp)
    end

    pure subroutine calculate_kernel_bandwidth(n, particles, weights, h)
        integer, intent(in) :: n
        real(dp), intent(in) :: particles(:,:), weights(:)
        real(dp), intent(out) :: h(size(particles,dim=1))
        integer :: i, ndims
        call debug_error_condition(size(particles,dim=2) /= n, 'mismatch in particles shape')
        call debug_error_condition(size(weights) /= n, 'number of weights must match number of particles')
        ndims = size(h)
        do concurrent (i=1:ndims)
            h(i) = std(particles(i,:), weights)
        end do
        h = max(1.0_dp, h)
    end

    impure subroutine apply_kernel_smoothing(h, h_scale, n, particles)
        real(dp), intent(in) :: h(:), h_scale
        integer, intent(in) :: n
        real(dp), intent(inout) :: particles(:,:)
        real(dp) :: new_particles(size(particles,dim=1),size(particles,dim=2)), cov_diag(size(h),size(h)), &
                    noise(size(particles,dim=1),size(particles,dim=2)), mu_zero(size(h))
        integer :: particles_ii(n), i
        call debug_error_condition(size(particles,dim=1) /= size(h), 'particles must be same dimension as h')
        call debug_error_condition(size(particles,dim=2) /= n, 'mismatch in particles shape')
        call random_uniform(particles_ii, 1, n)
        cov_diag = 0.0_dp
        do concurrent (i=1:size(h))
            cov_diag(i,i) = h(i)**2
        end do
        mu_zero = 0.0_dp
        call random_normal(noise, mu_zero, cov_diag)
        do concurrent (i=1:n)
            new_particles(:,i) = particles(:,particles_ii(i)) + h_scale*noise(:,i)
        end do
        particles = new_particles
    end

    pure subroutine fly_constant_acceleration(cart6_state, dt, max_spd)
        real(dp), intent(inout) :: cart6_state(:)
        real(dp), intent(in) :: dt, max_spd
        real(dp) :: spd_scale
        call debug_error_condition(size(cart6_state) /= 6, 'fly_constant_acceleration assumes state [x, y, vx, vy, ax, ay]')
        cart6_state(3:4) = cart6_state(3:4) + dt*cart6_state(5:6)
        spd_scale = max(1.0_dp, sqrt(cart6_state(3)**2 + cart6_state(4)**2)/max_spd)
        cart6_state(3:4) = cart6_state(3:4)/spd_scale
        cart6_state(1:2) = cart6_state(1:2) + dt*cart6_state(3:4)
        if (ballistic_particles .and. (cart6_state(2) <= 0.0_dp)) then
            cart6_state(2) = 0.0_dp
            cart6_state(4) = 0.0_dp
            cart6_state(6) = 0.0_dp
        end if
    end

    impure subroutine propagate_particles(dt, max_spd, max_acc, n, cartesian_particles, jerk_sig)
        real(dp), intent(in) :: dt, max_spd, max_acc
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:)
        real(dp), intent(in) :: jerk_sig(:)
        real(dp) :: jx(n), jy(n), acc_scale
        integer :: i
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(jerk_sig) /= n, 'mismatch in jerk_sig shape')
        if (ballistic_particles) then
            jx = 0.0_dp
            jy = 0.0_dp
        else
            call random_normal(jx, 0.0_dp, 1.0_dp)
            jx = jx*jerk_sig
            call random_normal(jy, 0.0_dp, 1.0_dp)
            jy = jy*jerk_sig
        end if
        do concurrent (i=1:n)
            if (ballistic_particles) then
                cartesian_particles(5,i) = 0.0_dp
                cartesian_particles(6,i) = -1.0_dp*g
            else
                cartesian_particles(5,i) = cartesian_particles(5,i) + dt*jx(i)
                cartesian_particles(6,i) = cartesian_particles(6,i) + dt*jy(i)
                acc_scale = max(1.0_dp, sqrt(cartesian_particles(5,i)**2 + cartesian_particles(6,i)**2)/max_acc)
                cartesian_particles(5:6,i) = cartesian_particles(5:6,i)/acc_scale
            end if
            call fly_constant_acceleration(cartesian_particles(:,i), dt, max_spd)
        end do
    end

  pure subroutine apply_constraints(obs_cart, max_rng, max_spd, max_acc, n, cartesian_particles, polar_particles)
        real(dp), intent(in) :: obs_cart(:), max_rng, max_spd, max_acc
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:)
        real(dp), intent(out) :: polar_particles(:,:)
        real(dp) :: los(2), ang, spd_scale, acc_scale
        integer :: i
        call debug_error_condition(size(cartesian_particles, dim=1) /= size(obs_cart), 'mismatch in obs_cart vs particle dims')
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(polar_particles, dim=2) /= n, 'mismatch in polar_particles shape')
        do concurrent (i=1:n)
            los = cartesian_particles(1:2,i) - obs_cart(1:2)
            if ((los(1)**2 + los(2)**2) > max_rng**2) then
                ang = atan2(los(2), los(1))
                cartesian_particles(1,i) = obs_cart(1) + max_rng*cos(ang)
                cartesian_particles(2,i) = obs_cart(2) + max_rng*sin(ang)
            end if
            spd_scale = max(1.0_dp, sqrt(cartesian_particles(3,i)**2 + cartesian_particles(4,i)**2)/max_spd)
            cartesian_particles(3:4,i) = cartesian_particles(3:4,i)/spd_scale
            acc_scale = max(1.0_dp, sqrt(cartesian_particles(5,i)**2 + cartesian_particles(6,i)**2)/max_acc)
            cartesian_particles(5:6,i) = cartesian_particles(5:6,i)/acc_scale
        end do
        call convert_particles_cart2pol(obs_cart, n, cartesian_particles, polar_particles)
    end

    pure function rmse(predicted, observed) result(val)
        real(dp), intent(in) :: predicted(:), observed
        real(dp) :: val
        real(dp) :: diff2(size(predicted))
        diff2 = (predicted - observed)**2
        val = sqrt(avg(diff2))
    end

    pure function neff(weights) result(val)
        real(dp), intent(in) :: weights(:)
        real(dp) :: val
        real(dp) :: weights2(size(weights))
        weights2 = weights**2
        val = 1.0_dp/dsum(weights2)
    end

    impure subroutine dump_particles(filename, t, ref_cart, cartesian_particles, weights)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: t, ref_cart(:), cartesian_particles(:,:), weights(:)
        real(dp) :: est_cart(size(ref_cart))
        integer :: fid, i, n
        character(len=1024) :: fmtstr
        call debug_error_condition(size(cartesian_particles,dim=1) /= size(ref_cart), &
                                   'state estimate dimensions need to match particle dimensions')
        call debug_error_condition(size(cartesian_particles,dim=2) /= size(weights), &
                                   'number of particles needs to match number of weights')
        n = size(cartesian_particles, dim=2)
        call generate_state_estimate(est_cart, n, cartesian_particles, weights)
        open(newunit=fid, file=filename, action='write')
        write(fid,'(a)') 'number,time_sec,x_ft,y_ft,vx_fps,vy_fps,ax_fps2,ay_fps2,weight'
        write(fmtstr,'(a,i0,a)') '(i0,",",f0.1,",",',size(est_cart),'(e13.6,","),a)'
        write(unit=fid, fmt=trim(fmtstr)) 0, -1.0, ref_cart, 'reference'
        write(unit=fid, fmt=trim(fmtstr)) 0, -1.0, est_cart, 'estimate'
        write(fmtstr,'(a,i0,a)') '(i0,",",f0.1,",",',size(est_cart),'(e13.6,","),e13.6)'
        do i=1,n
            write(unit=fid, fmt=trim(fmtstr)) i, t, cartesian_particles(:,i), weights(i)
        end do
        close(fid)
    end

end module pf


program ex_particle_filter
use, non_intrinsic :: pf
implicit none

    integer, parameter :: max_trials = 1024
    real(dp), parameter :: max_rng   = 500.0_dp*nmi2ft, &
                           max_spd   = 10000.0_dp, &
                           max_acc   = 9.0_dp*g, &
                           dt        = 1.0_dp, &
                           meas_sig(3) = [10.0_dp*nmi2ft, 5.0_dp*deg2rad_dp, 200.0_dp] ! poor measurements
!                           meas_sig(3) = [100.0_dp, 0.1_dp*deg2rad_dp, 10.0_dp] ! standard measurements
!                           meas_sig(3) = [1.0_dp, 0.001_dp, 1.0_dp] ! exquisite measurements
!                           meas_sig(3) = [-1.0_dp, 0.1_dp*deg2rad_dp, -1.0_dp] ! bearing only

    procedure(resample_cartesian_particles), pointer :: resample_subroutine
    procedure(calculate_weights), pointer :: calculate_weights_subroutine
    integer :: dr, tmax_int, spd, rough_fac_int, neff_pass_int, num_particles, trial, num_particles_ii, init_meas_sig, &
               resampling_method, jerk_sig_lo_int, jerk_sig_hi_int, roughening_method, fid, init_vt_method, likelihood_method
    real(dp) :: obs_cart(6), tgt_cart(6), tgt_pol(3), t, tmax, rough_fac, neff_thresh, neff_pass, neff0, &
                meas(3), est_cart(6), est_pol(3), meas_sig_fac, jerk_sig_lo, jerk_sig_hi, h(size(est_cart)), &
                x_err(max_trials), y_err(max_trials), pos_err(max_trials), &
                vx_err(max_trials), vy_err(max_trials), vel_err(max_trials), &
                ax_err(max_trials), ay_err(max_trials), acc_err(max_trials), &
                rng_err(max_trials), ang_err(max_trials), rngrt_err(max_trials)
    real(dp), allocatable :: cartesian_particles(:,:), weights(:), polar_particles(:,:), jerk_sig(:)
    character(len=1024) :: fmtstr, resample_strategy, run_description, num_particles_str, init_meas_sig_str, &
                           neff_pass_int_str, rough_fac_int_str, resampling_method_str, tmax_int_str, dr_str, spd_str, &
                           max_trials_str, trial_str, filename, jerk_sig_lo_int_str, jerk_sig_hi_int_str, roughening_method_str, &
                           roughen_strategy, init_vt_mode, init_vt_method_str, likelihood_method_str, likelihood_mode

    open(newunit=fid, file='/valinor/last-run.csv', action='write')
    write(fid,'(a)') 'max_trials,num_particles,init_meas_sig,vt_mode,likelihood_method,resampling_method,roughening_method,'// &
                     'rough_fac,neff_thresh,'// &
                     'jerk_sig_lo,jerk_sig_hi,tmax_sec,downrange_nmi,spd_mach,rmse_rng,rmse_ang,rmse_rngrt,rmse_x,rmse_y,'// &
                     'rmse_vx,rmse_vy,rmse_ax,rmse_ay,rmse_pos,rmse_vel,rmse_acc,rmse_composite,max_pos_err'
    write(max_trials_str,'(a,i0)') 'max_trials: ',max_trials
    do num_particles_ii=11,11,10
        num_particles = 10**(num_particles_ii/10)*1000*mod(num_particles_ii,10)
        if (num_particles == 0) cycle
        if (allocated(cartesian_particles)) deallocate(cartesian_particles)
        if (allocated(weights)) deallocate(weights)
        if (allocated(polar_particles)) deallocate(polar_particles)
        if (allocated(jerk_sig)) deallocate(jerk_sig)
        allocate(cartesian_particles(size(obs_cart),num_particles), &
                 weights(num_particles), &
                 polar_particles(3,num_particles), &
                 jerk_sig(num_particles))
        write(num_particles_str,'(a,i0)') ' | num_particles: ',num_particles
    do neff_pass_int=5,100,5 !! respample when Neff is 5-100% of num_particles                                        
        neff_thresh = real(neff_pass_int, kind=dp)/100.0_dp                                                           
        neff_pass = ceiling(neff_thresh*num_particles)
        write(neff_pass_int_str,'(a,f0.2)') ' | neff_thresh: ',neff_thresh
    do rough_fac_int=5,100,5 !! scaling Neff**(-1/(d+4))
        rough_fac = real(rough_fac_int,dp)/100.0_dp
        write(rough_fac_int_str,'(a,f0.2)') ' | rough_fac: ',rough_fac
    do init_vt_method=1,1 !! tangential velocity component initialization
        select case (init_vt_method)
            case (1) !! Vtotal uniform from Vradial to Vmax
                vt_initialization_mode = 1
                init_vt_mode = 'uniform'
            case (2) !! Vtotal log-uniform from Vradial to Vmax
                vt_initialization_mode = 2
                init_vt_mode = 'log_uniform'
            case (3) !! Vtotal normal about Vradial with sigma (Vmax-Vradial)/1
                vt_initialization_mode = 3
                init_vt_mode = 'normal_sig_1'
            case (4) !! Vtotal normal about Vradial with sigma (Vmax-Vradial)/2
                vt_initialization_mode = 4
                init_vt_mode = 'normal_sig_2'
            case (5) !! Vtotal normal about Vradial with sigma (Vmax-Vradial)/3
                vt_initialization_mode = 5
                init_vt_mode = 'normal_sig_3'
            case default
                error stop 'not implemented'
        end select
        write(init_vt_method_str,'(a)') ' | vt_method: '//trim(init_vt_mode)
    do likelihood_method=1,2 !! gaussian or student_t
        select case (likelihood_method)
            case (1)
                likelihood_mode = 'gaussian'
                calculate_weights_subroutine => calculate_weights_gaussian
            case (2)
                likelihood_mode = 'student_t'
                calculate_weights_subroutine => calculate_weights_student_t
            case default
                error stop 'not implemented'
        end select
        write(likelihood_method_str,'(a)') ' | likelihood_method: '//trim(likelihood_mode)
    do resampling_method=0,0 !! select case on method used
        select case (resampling_method)
            case (0)
                resample_strategy = 'covariance_regularization'
            case (1)
                resample_strategy = 'systematic'
                resample_subroutine => resample_systematic
            case (2)
                resample_strategy = 'residual'
                resample_subroutine => resample_residual
            case (3)
                resample_strategy = 'stratified'
                resample_subroutine => resample_stratified
            case (4)
                resample_strategy = 'multinomial'
                resample_subroutine => resample_multinomial
            case default
                error stop 'not implemented yet'
        end select
        write(resampling_method_str,'(a)') ' | resample_strategy: '//trim(resample_strategy)
    do roughening_method=1,2 !! select case on method used
        select case (roughening_method)
            case (1)
                roughen_strategy = 'cov_reg'
            case (2)
                roughen_strategy = 'neff_cov_reg'
            case (3)
                roughen_strategy = 'kern_smooth'
            case (4)
                roughen_strategy = 'neff_kern_smooth'
            case (5)
                roughen_strategy = 'jitter_vxy'
            case (6)
                roughen_strategy = 'neff_jitter_vxy'
            case (7)
                roughen_strategy = 'jitter_axy'
            case (8)
                roughen_strategy = 'neff_jitter_axy'
            case default
                error stop 'not implemented yet, want to add vxy jitter and axy jitter options'
        end select
        write(roughening_method_str,'(a)') ' | roughen_strategy: '//trim(roughen_strategy)
    do init_meas_sig=1,1,-1 !! multiply meas_sig by init_meas_sig for particle initialization
        meas_sig_fac = 1.0_dp/real(init_meas_sig, kind=dp)
        write(init_meas_sig_str,'(a,f0.4)') ' | init_meas_sig: ',meas_sig_fac
    do jerk_sig_lo_int=-325,-325 !! -325,-325 ! no args defaults
        jerk_sig_lo = 10.0_dp**(real(jerk_sig_lo_int,dp)/100.0_dp)
        write(jerk_sig_lo_int_str,'(a,e13.6)') ' | jerk_sig_lo: ',jerk_sig_lo
    do jerk_sig_hi_int=-75,-75 !! -75,-75 ! no args defaults
        jerk_sig_hi = 10.0_dp**(real(jerk_sig_hi_int,dp)/100.0_dp)
        write(jerk_sig_hi_int_str,'(a,e13.6)') ' | jerk_sig_hi: ',jerk_sig_hi
    do tmax_int=201,201 !! performance should increase with tmax until target passes max_spd
        tmax = real(tmax_int, kind=dp)
        write(tmax_int_str,'(a,f0.1)') ' | tmax: ',tmax
        if (mod(tmax_int, 2) == 1) then
            ballistic_particles = .true. !! odd tmax assume ballistic ay
        else
            ballistic_particles = .false. !! even tmax assume random jerk
        end if
    do dr=50,250,50 !! NMI
        write(dr_str,'(a,i0)') ' | dr: ',dr
    loop_spd: do spd=1500,9000,1500 !! ~Mach
        write(spd_str,'(a,i0)') ' | spd: ',spd
    !$omp parallel do default(firstprivate) private(neff0, t) shared(x_err, y_err, vx_err, vy_err, ax_err, ay_err, pos_err, &
    !$omp&                                                           vel_err, acc_err, rng_err, ang_err, rngrt_err)
    loop_trials: do trial=1,max_trials
        !! reset starting positions for each run
        obs_cart = 0.0_dp
        tgt_cart = 0.0_dp !! initialize all state components to zero (0.0)
        tgt_cart(1) = dr*nmi2ft                                       !! x position [ft]
        tgt_cart(2) = 1.0_dp                                          !! y position [ft]
        tgt_cart(3) = -real(spd,dp)*cos(45.0_dp*deg2rad_dp) !! vx velocity [ft/sec]
        tgt_cart(4) = real(spd,dp)*sin(45.0_dp*deg2rad_dp)  !! vy velocity [ft/sec]
        tgt_cart(5) = 0.0_dp                                          !! ax acceleration [ft/sec**2]
        tgt_cart(6) = -1.0_dp*g                                       !! ay acceleration [ft/sec**2]
        call perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
        call generate_measurements(obs_cart, tgt_cart, meas_sig, max_rng, max_spd, 1, meas)
        call initialize_particles(obs_cart, meas, meas_sig_fac*meas_sig, max_rng, max_spd, &
                                  num_particles, cartesian_particles, weights, jerk_sig)
        if (debug_init) then
            write(filename,'(i0,a,i0,a)') num_particles,'-particles-',init_meas_sig,'-sig-init.csv'
            call dump_particles(trim(filename), t, tgt_cart, cartesian_particles, weights)
            cycle
        end if
        write(trial_str,'(a,i0)') ' | trial: ',trial
        run_description = trim(max_trials_str)//trim(trial_str)//trim(num_particles_str)//trim(init_vt_method_str)// &
                          trim(likelihood_method_str)//trim(init_meas_sig_str)// &
                          trim(neff_pass_int_str)//trim(jerk_sig_lo_int_str)//trim(jerk_sig_hi_int_str)// &
                          trim(resampling_method_str)//trim(roughening_method_str)//trim(rough_fac_int_str)//trim(tmax_int_str)// &
                          trim(dr_str)//trim(spd_str)
        if (debug) then
            write(*,'(a)') trim(run_description)
            write(*,'(6(a,f0.1))') 'TARGET TRUTH x: ',tgt_cart(1),', y: ',tgt_cart(2),', vx: ',tgt_cart(3), &
                                             ', vy: ',tgt_cart(4), ', ax: ',tgt_cart(5),', ay: ',tgt_cart(6)
            call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
            write(*,'(a,f0.1)') '  Neff at initialization: ',neff(weights)
            write(*,'(6(a,f0.1))') '    ESTIMATE x: ',est_cart(1),', y: ',est_cart(2),', vx: ',est_cart(3), &
                                             ', vy: ',est_cart(4), ', ax: ',est_cart(5),', ay: ',est_cart(6)
            write(*,'(a)') ''
        end if
        t = 0.0_dp
        do while ((tgt_cart(2) > 0.0_dp) .and. (t < tmax))
            !! advance time and target independent of anything else
            t = t + dt
            call fly_constant_acceleration(tgt_cart, dt, huge(max_spd))
            call fly_constant_acceleration(obs_cart, dt, huge(max_spd))
            if (debug) then
                write(*,'(a,f0.1)') 't=',t
                write(*,'(6(a,f0.1))') 'TARGET TRUTH x: ',tgt_cart(1),', y: ',tgt_cart(2),', vx: ',tgt_cart(3), &
                                                 ', vy: ',tgt_cart(4), ', ax: ',tgt_cart(5),', ay: ',tgt_cart(6)
            end if
            !! update tgt_pol
            call perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
            !! generate new measurement from observer perspective
            call generate_measurements(obs_cart, tgt_cart, meas_sig, max_rng, max_spd, 1, meas)
            !! propagate current particles
            call propagate_particles(dt, max_spd, max_acc, num_particles, cartesian_particles, jerk_sig)
            if (debug) then
                call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
                write(*,'(6(a,f0.1))') '   POST-PROP x: ',est_cart(1),', y: ',est_cart(2),', vx: ',est_cart(3), &
                                                 ', vy: ',est_cart(4), ', ax: ',est_cart(5),', ay: ',est_cart(6)
            end if
            !! apply constraints and convert particles from cartesian to polar representation
            call apply_constraints(obs_cart, max_rng, max_spd, max_acc, num_particles, cartesian_particles, polar_particles)
            if (debug) then
                call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
                write(*,'(6(a,f0.1))') ' CONSTRAINED x: ',est_cart(1),', y: ',est_cart(2),', vx: ',est_cart(3), &
                                                 ', vy: ',est_cart(4), ', ax: ',est_cart(5),', ay: ',est_cart(6)
            end if
            !! calculate weights, resample if Neff is below neff_pass (acceptable percentage of original num_particles)
            call calculate_weights_subroutine(meas, meas_sig, num_particles, polar_particles, weights)
            neff0 = neff(weights)
            if (debug) then
                write(*,'(a,f0.1)') '  Neff before resample: ',neff0
            end if
            if (neff0 < neff_pass) then
!                call resample_subroutine(num_particles, cartesian_particles, weights, jerk_sig)
!                if (debug) then
!                    call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
!                    write(*,'(6(a,f0.1))') 'POSTRESAMPLE x: ',est_cart(1),', y: ',est_cart(2),', vx: ',est_cart(3), &
!                                                     ', vy: ',est_cart(4), ', ax: ',est_cart(5),', ay: ',est_cart(6)
!                end if
!                if (roughening_method > 2) call calculate_kernel_bandwidth(num_particles, cartesian_particles, weights, h)
                select case (roughening_method)
                    case (1) !! covariance regularization
                        call regularize_covariance(rough_fac, num_particles, cartesian_particles, weights, jerk_sig)
                    case (2) !! covariance regularization with neff scaling
                        call regularize_covariance(rough_fac*neff0**(-1.0_dp/(size(cartesian_particles,dim=1)+4)), &
                                              num_particles, cartesian_particles, weights, jerk_sig)
                    case (3) !! pre-resample bandwidth kernel smoothing
                        call apply_kernel_smoothing(h, rough_fac, num_particles, cartesian_particles)
                    case (4) !! pre-resample bandwidth kernel smoothing with neff scaling
                        call apply_kernel_smoothing(h, rough_fac*neff0**(-1.0_dp/(size(cartesian_particles,dim=1)+4)), &
                                                    num_particles, cartesian_particles)
                    case (5) !! pre-resample bandwidth kernel smoothing only vx and vy
                        h(1) = 0.0_dp; h(2) = 0.0_dp; h(5) = 0.0_dp; h(6) = 0.0_dp;
                        call apply_kernel_smoothing(h, rough_fac, num_particles, cartesian_particles)
                    case (6) !! pre-resample bandwidth kernel smoothing only vx and vy with neff scaling
                        h(1) = 0.0_dp; h(2) = 0.0_dp; h(5) = 0.0_dp; h(6) = 0.0_dp;
                        call apply_kernel_smoothing(h, rough_fac*neff0**(-1.0_dp/(size(cartesian_particles,dim=1)+4)), &
                                                    num_particles, cartesian_particles)
                    case (7) !! pre-resample bandwidth kernel smoothing only ax and ay
                        h(1) = 0.0_dp; h(2) = 0.0_dp; h(3) = 0.0_dp; h(4) = 0.0_dp;
                        call apply_kernel_smoothing(h, rough_fac, num_particles, cartesian_particles)
                    case (8) !! pre-resample bandwidth kernel smoothing only ax and ay with neff scaling
                        h(1) = 0.0_dp; h(2) = 0.0_dp; h(3) = 0.0_dp; h(4) = 0.0_dp;
                        call apply_kernel_smoothing(h, rough_fac*neff0**(-1.0_dp/(size(cartesian_particles,dim=1)+4)), &
                                                    num_particles, cartesian_particles)
                    case default
                        error stop 'none implemented'
                end select
                call resample_jerk(jerk_sig, jerk_sig_lo, jerk_sig_hi)
                if (debug) then
                    call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
                    write(*,'(6(a,f0.1))') '  POST-ROUGH x: ',est_cart(1),', y: ',est_cart(2),', vx: ',est_cart(3), &
                                                     ', vy: ',est_cart(4), ', ax: ',est_cart(5),', ay: ',est_cart(6)
                end if
            end if
            if (debug) then
                write(*,'(a)') ''
                if (.not.all(nearly(est_cart, est_cart))) error stop 'nan detected'
            end if
        end do
        call generate_state_estimate(est_pol, num_particles, polar_particles, weights)
        rng_err(trial) = tgt_pol(1) - est_pol(1)
        ang_err(trial) = mod(tgt_pol(2) - est_pol(2) + pi_dp, twopi_dp) - pi_dp
        rngrt_err(trial) = tgt_pol(3) - est_pol(3)
        call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
        x_err(trial) = tgt_cart(1) - est_cart(1)
        y_err(trial) = tgt_cart(2) - est_cart(2)
        vx_err(trial) = tgt_cart(3) - est_cart(3)
        vy_err(trial) = tgt_cart(4) - est_cart(4)
        ax_err(trial) = tgt_cart(5) - est_cart(5)
        ay_err(trial) = tgt_cart(6) - est_cart(6)
        pos_err(trial) = sqrt((tgt_cart(1) - est_cart(1))**2 + (tgt_cart(2) - est_cart(2))**2)
        vel_err(trial) = sqrt((tgt_cart(3) - est_cart(3))**2 + (tgt_cart(4) - est_cart(4))**2)
        acc_err(trial) = sqrt((tgt_cart(5) - est_cart(5))**2 + (tgt_cart(6) - est_cart(6))**2)
        if (debug_run) then
            fmtstr = '(4(a,e13.6))'
            write(*,fmtstr) trim(run_description)//' | pos_err: ',pos_err(trial),' | vel_err: ',vel_err(trial), &
                ' | acc_err: ',acc_err(trial),' | COMPOSITE: ',geomean([pos_err(trial), vel_err(trial), acc_err(trial)])
        end if
    end do loop_trials
    !$omp end parallel do
        if (debug_init) cycle
        !! trials done
        fmtstr = '(2(i0,","),e13.6,",",4(a,","),2(f0.4,","),2(e13.6,","),f0.1,",",2(i0,","),13(e13.6,","),e13.6)'
        if (all(nearly(pos_err, pos_err))) then
            write(fid,fmt=trim(fmtstr)) max_trials,num_particles,meas_sig_fac,trim(init_vt_mode),trim(likelihood_mode), &
                                        trim(resample_strategy), &
                                        trim(roughen_strategy),rough_fac,neff_thresh,jerk_sig_lo,jerk_sig_hi,tmax,dr,spd, &
                                        rmse(rng_err,0.0_dp),rmse(ang_err,0.0_dp),rmse(rngrt_err,0.0_dp),rmse(x_err,0.0_dp), &
                                        rmse(y_err,0.0_dp),rmse(vx_err,0.0_dp),rmse(vy_err,0.0_dp),rmse(ax_err,0.0_dp), &
                                        rmse(ay_err,0.0_dp),rmse(pos_err,0.0_dp),rmse(vel_err,0.0_dp),rmse(acc_err,0.0_dp), &
                                        geomean([rmse(pos_err,0.0_dp),rmse(vel_err,0.0_dp),rmse(acc_err,0.0_dp)]),maxval(pos_err)
            flush(fid)
        else
            error stop 'nan should not happen'
        end if
        fmtstr = '(a,e13.6)'
        run_description = trim(max_trials_str)//trim(num_particles_str)//trim(init_meas_sig_str)//trim(neff_pass_int_str)// &
                          trim(rough_fac_int_str)//trim(init_vt_method_str)//trim(likelihood_method_str)// &
                          trim(resampling_method_str)// &
                          trim(roughening_method_str)
        write(*,fmtstr) trim(run_description)//' | RMSE position: ',rmse(pos_err,0.0_dp)
        if (debug) then
        call sort(pos_err)
            do trial=1,9
                write(*,'(i0,a,f0.1,a)') trial*10,'th percentile pos_err: ', &
                                         pos_err(ceiling(real(trial,dp)/10.0_dp*max_trials)),' ft'
            end do
        end if
    end do loop_spd
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do

    close(fid)

end program ex_particle_filter
