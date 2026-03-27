module eoir_uav_mod
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: nmi2ft_dp, ft2nmi_dp, deg2rad_dp, km2nmi_dp
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: math_utils, only: linspace
use, non_intrinsic :: statistics, only: normalize_min_max
implicit none
private

    real(dp), parameter :: R_air = 1716.46_dp !! https://en.wikipedia.org/wiki/Gas_constant#Specific_gas_constant

    public :: dp, km2nmi_dp, endurance_hr, cost_m, leg_hr, calc_patches, linspace, normalize_min_max

contains

    pure elemental function endurance_hr(mach, alt_kft) result(val)
        real(dp), intent(in) :: mach, alt_kft
        real(dp) :: val
        val = -18.75_dp*mach**2 + 8.0893_dp*mach + 0.01_dp*alt_kft**2 + 0.05_dp*alt_kft + 9.2105_dp
    end function endurance_hr

    pure elemental function cost_m(mach, alt_kft) result(val)
        real(dp), intent(in) :: mach, alt_kft
        real(dp) :: val
        val = 50.0_dp*mach**2 - 35.0_dp*mach + 0.03_dp*alt_kft**2 - 0.2_dp*alt_kft + 11.0_dp
    end function cost_m

    pure elemental function imperfect_gamma(T_R) result(val) !! https://www.grc.nasa.gov/www/BGH/realspec.html
        real(dp), intent(in) :: T_R
        real(dp) :: val
        real(dp) :: Theta_T_R
        call debug_error_condition(T_R <= 0.0_dp, 'module EOIR_UAV_MOD :: funciton imperfect_gamma invalid for T_R <= 0.0')
        Theta_T_R = 5500.0_dp/T_R
        val = 1.0_dp + 0.4_dp/(1.0_dp + 0.4_dp*(Theta_T_R**2*exp(Theta_T_R)/(exp(Theta_T_R)-1.0_dp)**2))
    end function imperfect_gamma

    pure elemental function t_alt(alt_ft) result(val) !! https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/atmos.html
        real(dp), intent(in) :: alt_ft
        real(dp) :: val
        if (alt_ft < 36152.0_dp) then
            val = 59.0_dp - 0.00356_dp*alt_ft
        else if(alt_ft > 82345.0_dp) then
            val = -205.05_dp + 0.00164_dp*alt_ft
        else !! 36152.0 <= alt_ft <= 82345
            val = -70.0_dp
        end if
        val = val + 459.67_dp
    end function t_alt

    pure elemental function mach1(alt_ft) result(val)
        real(dp), intent(in) :: alt_ft
        real(dp) :: val
        real(dp) :: T_R_at_alt, y_at_alt
        T_R_at_alt = t_alt(alt_ft)
        y_at_alt = imperfect_gamma(T_R_at_alt)
        val = sqrt(y_at_alt*R_air*T_R_at_alt) ! ft/sec
    end function mach1

    pure elemental function leg_hr(leg_nmi, mach, alt_kft) result(val)
        real(dp), intent(in) :: leg_nmi, mach, alt_kft
        real(dp) :: val
        val = leg_nmi*nmi2ft_dp/(mach1(alt_kft*1000.0_dp)*mach)/3600.0_dp
    end function leg_hr

    pure elemental subroutine calc_patches(mach, alt_kft, fov_az_hw_deg, fov_el_hw_deg, az_center_deg, el_center_deg, &
                                           mission_width_nmi, mission_length_nmi, platform_mission_time_hr, &
                                           patch_area_nmi2, min_platforms_required, max_total_survey_area_nmi2, &
                                           platform_route_nmi, platform_route_hr)
        real(dp), intent(in) :: mach, alt_kft, fov_az_hw_deg, fov_el_hw_deg, az_center_deg, el_center_deg, &
                                mission_width_nmi, mission_length_nmi, platform_mission_time_hr
        real(dp), intent(out) :: patch_area_nmi2
        integer, intent(out) :: min_platforms_required
        real(dp), intent(out) :: max_total_survey_area_nmi2, platform_route_nmi, platform_route_hr
        real(dp) :: patch_width_nmi, patch_length_nmi, min_strip_length_nmi, min_strip_time_hr, strip_area_nmi2, mission_area_nmi2
        integer :: nrow, min_strips_required
        call debug_error_condition(.not.(nearly(az_center_deg, 0.0_dp) .and. nearly(el_center_deg, 0.0_dp)), &
                                   'module EOIR_UAV_MOD :: function patch_area_nmi2 only implemented for pure lookdown case')
        patch_length_nmi = alt_kft*1000.0_dp*ft2nmi_dp*tan(fov_el_hw_deg*deg2rad_dp)*2.0_dp
        patch_width_nmi = alt_kft*1000.0_dp*ft2nmi_dp*tan(fov_az_hw_deg*deg2rad_dp)*2.0_dp
        patch_area_nmi2 = patch_length_nmi*patch_width_nmi
        nrow = ceiling(mission_length_nmi/patch_length_nmi) ! round up for 100% coverage in length
        !! platform must return to mission area starting side
        !! minimum mission would be: entry to mission area 0.5*patch_length_nmi
        !!                           fly nrow-1 patches (n-1)*patch_length_nmi
        !!                           fly laterally patch_width_nmi
        !!                           fly nrow-1 patches (n-1)*patch_length_nmi
        !!                           exit mission area 0.5*patch_length_nmi
        min_strip_length_nmi = 0.5_dp*patch_length_nmi + & ! entry to mission area
                               (nrow-1)*patch_length_nmi + patch_width_nmi + (nrow-1)*patch_length_nmi + & ! down mission_length_nmi, over patch_width_nmi, and back mission_length_nmi
                               0.5_dp*patch_length_nmi ! exit mission area
        min_strip_time_hr = leg_hr(min_strip_length_nmi, mach, alt_kft)
        if (min_strip_time_hr < platform_mission_time_hr) then !! valid configuration, proceed to determine min_platforms_required, platform_route_nmi, and platform_route_hr
            strip_area_nmi2 = patch_area_nmi2*(2*nrow) ! each strip returns platform to mission entry side
            min_platforms_required = 0
            platform_route_hr = huge(1.0_dp)
            do while (platform_route_hr > platform_mission_time_hr)
                min_platforms_required = min_platforms_required + 1
                mission_area_nmi2 = mission_length_nmi*mission_width_nmi/min_platforms_required ! area per platform = length * width / platforms
                min_strips_required = ceiling(mission_area_nmi2/strip_area_nmi2) ! round up for 100% coverage in width
                platform_route_nmi = patch_length_nmi + & ! entry and exit from mission area
                                     (2*(nrow - 1)*patch_length_nmi + patch_width_nmi)*min_strips_required + & ! individual strip length
                                     patch_width_nmi*(min_strips_required - 1) ! movement between strips
                platform_route_hr = leg_hr(platform_route_nmi, mach, alt_kft)
            end do
            !! maximum surveyed area = patch_width_nmi * distance traveled in platform_mission_time_hr * num_platforms
            max_total_survey_area_nmi2 = patch_width_nmi*platform_mission_time_hr*3600.0_dp* &
                                            (mach1(alt_kft*1000.0_dp)*mach*ft2nmi_dp)*min_platforms_required
        else !! invalid configuration, set primary metrics to -1 and exit
            call debug_error_condition(.true., 'invalid configuration scoring not implemented')
            min_platforms_required = -1
            platform_route_nmi = -1.0_dp
            platform_route_hr = -1.0_dp
            max_total_survey_area_nmi2 = -1.0_dp
        end if
    end subroutine calc_patches

end module eoir_uav_mod


program main
use, non_intrinsic :: eoir_uav_mod
implicit none

    integer, parameter :: mach_list_size = 11, alt_list_size = 11
    real(dp), parameter :: ingress_nmi = 100.0_dp*km2nmi_dp, egress_nmi = 100.0_dp*km2nmi_dp, &
                           mission_width_nmi = 100.0_dp*km2nmi_dp, mission_length_nmi = 100.0_dp*km2nmi_dp, &
                           fov_deg_list(*) = [15.0_dp, 30.0_dp, 60.0_dp], sensor_cost_m_list(*) = [0.05_dp, 1.0_dp, 10.0_dp], &
                           min_mach = 0.40_dp, max_mach = 0.90_dp, min_alt_kft = 5.0_dp, max_alt_kft = 25.0_dp

    integer :: mach_ii, alt_ii, fov_ii, fid, c_ii, min_platforms_required
    real(dp) :: mach, alt_kft, fov_deg, platform_endurance_hr, platform_ingress_time_hr, platform_egress_time_hr, &
                platform_mission_time_hr, &
                platform_mission_patch_nmi2, airframe_cost_m, sensor_cost_m, platform_cost_m, &
                platform_mission_route_nmi, platform_mission_route_hr, max_total_survey_area_nmi2, &
                max_total_survey_area_nmi2_per_cost_m, &
                mission_reps, mission_reps_per_cost_m, mission_reps_per_hour, mission_reps_per_hour_per_cost_m, &
                mission_one_rep_time_hr, cost_m_per_mission_one_rep_time_hr, total_cost_m, &
                mach_list(mach_list_size), alt_list(alt_list_size), results(37,mach_list_size*alt_list_size*size(fov_deg_list)), &
                absolute_performance_weight, efficiency_weight, total_weight
    character(len=128) :: fmt_str

    if (command_argument_count() == 2) then
        call get_command_argument(1, fmt_str)
        read(fmt_str,*) absolute_performance_weight
        call get_command_argument(2, fmt_str)
        read(fmt_str,*) efficiency_weight
    else
        absolute_performance_weight = 0.50_dp
        efficiency_weight = 0.50_dp
    end if
    total_weight = absolute_performance_weight + efficiency_weight

    call linspace(mach_list, min_mach, max_mach)
    call linspace(alt_list, min_alt_kft, max_alt_kft)

    !$omp parallel do default(firstprivate) shared(results)
    do mach_ii=1,mach_list_size
        mach = mach_list(mach_ii)
        do alt_ii=1,alt_list_size
            alt_kft = alt_list(alt_ii)
            platform_endurance_hr = endurance_hr(mach, alt_kft)
            platform_ingress_time_hr = leg_hr(ingress_nmi, mach, alt_kft)
            platform_egress_time_hr = leg_hr(egress_nmi, mach, alt_kft)
            platform_mission_time_hr = platform_endurance_hr - platform_ingress_time_hr - platform_egress_time_hr
            airframe_cost_m = cost_m(mach, alt_kft)
            do fov_ii=1,size(fov_deg_list)
                c_ii = fov_ii + (alt_ii - 1)*size(fov_deg_list) + (mach_ii - 1)*alt_list_size*size(fov_deg_list)
                fov_deg = fov_deg_list(fov_ii)
                sensor_cost_m = sensor_cost_m_list(fov_ii)
                platform_cost_m = airframe_cost_m + sensor_cost_m
                call calc_patches(mach, alt_kft, fov_deg*0.5_dp, fov_deg*0.5_dp, 0.0_dp, 0.0_dp, &
                                  mission_width_nmi, mission_length_nmi, platform_mission_time_hr, &
                                  platform_mission_patch_nmi2, min_platforms_required, max_total_survey_area_nmi2, &
                                  platform_mission_route_nmi, platform_mission_route_hr)
                total_cost_m = platform_cost_m*min_platforms_required
                mission_reps = platform_mission_time_hr/platform_mission_route_hr
                mission_reps_per_cost_m = mission_reps/total_cost_m
                mission_reps_per_hour = 1.0_dp/platform_mission_route_hr
                mission_reps_per_hour_per_cost_m = mission_reps_per_hour/total_cost_m
                mission_one_rep_time_hr = platform_ingress_time_hr + platform_mission_route_hr + platform_egress_time_hr
                cost_m_per_mission_one_rep_time_hr = total_cost_m/mission_one_rep_time_hr
                max_total_survey_area_nmi2_per_cost_m = max_total_survey_area_nmi2/total_cost_m
                !! fill results
                results( 1,c_ii) = mach
                results( 2,c_ii) = alt_kft
                results( 3,c_ii) = fov_deg
                results( 4,c_ii) = platform_ingress_time_hr
                results( 5,c_ii) = platform_egress_time_hr
                results( 6,c_ii) = platform_mission_time_hr
                results( 7,c_ii) = platform_endurance_hr
                results( 8,c_ii) = platform_mission_patch_nmi2
                results( 9,c_ii) = platform_mission_route_nmi
                results(10,c_ii) = platform_mission_route_hr
                results(11,c_ii) = airframe_cost_m
                results(12,c_ii) = sensor_cost_m
                results(13,c_ii) = platform_cost_m
                results(14,c_ii) = real(min_platforms_required, kind=dp)
                results(15,c_ii) = total_cost_m
                !! results(16,c_ii) = norm_total_cost_m
                results(17,c_ii) = mission_reps
                results(18,c_ii) = mission_reps_per_cost_m
                results(19,c_ii) = mission_reps_per_hour
                results(20,c_ii) = mission_reps_per_hour_per_cost_m
                results(21,c_ii) = mission_one_rep_time_hr
                results(22,c_ii) = cost_m_per_mission_one_rep_time_hr
                results(23,c_ii) = max_total_survey_area_nmi2
                results(24,c_ii) = max_total_survey_area_nmi2_per_cost_m
                results(36,c_ii) = sqrt(total_cost_m**2 + mission_one_rep_time_hr**2)
            end do
        end do
    end do
    !$omp end parallel do

    call normalize_min_max(results(15,:), results(16,:)) !! lower = better, total_cost_m

    call normalize_min_max(results(17,:), results(25,:)) !! higher = better, mission_reps
    call normalize_min_max(results(18,:), results(26,:)) !! higher = better, mission_reps_per_cost_m
    call normalize_min_max(results(19,:), results(27,:)) !! higher = better, mission_reps_per_hour
    call normalize_min_max(results(20,:), results(28,:)) !! higher = better, mission_reps_per_hour_per_cost_m

    call normalize_min_max(results(21,:), results(29,:)) !! lower = better, mission_one_rep_time_hr
    results(37,:) = sqrt(results(16,:)**2 + results(29,:)**2) !! RSS normalized total_cost_m and normalized mission_one_rep_time_hr to compare differing scales of values (cost ($M) and time (hr))

    results(29,:) = 1.0_dp - results(29,:) !! flip 0-1 so lower=better can be averaged with other higher=better scores

    call normalize_min_max(results(22,:), results(30,:)) !! lower = better, cost_m_per_mission_one_rep_time_hr
    results(30,:) = 1.0_dp - results(30,:) !! flip 0-1 so lower=better can be averaged with other higher=better scores

    call normalize_min_max(results(23,:), results(31,:)) !! higher = better, max_total_survey_area_nmi2
    call normalize_min_max(results(24,:), results(32,:)) !! higher = better, max_total_survey_area_nmi2_per_cost_m

    results(33,:) = (results(25,:) + results(27,:) + results(29,:) + results(31,:))/4.0_dp !! average absolute scores
    results(34,:) = (results(26,:) + results(28,:) + results(30,:) + results(32,:))/4.0_dp !! average efficiency scores
    results(35,:) = (absolute_performance_weight*results(33,:) + efficiency_weight*results(34,:))/total_weight !! weighted average total scores

    open(newunit=fid, file='/valinor/eoir-uav-trade-study.csv', action='write')
    write(fid,'(a,f0.4,a,f0.4,a)') 'mach,alt_kft,sensor_fov_deg,'// &
                     'platform_ingress_time_hr,platform_egress_time_hr,platform_mission_time_hr,platform_endurance_hr,'// &
                     'platform_survey_patch_nmi2,platform_mission_route_nmi,platform_mission_route_hr,'// &
                     'airframe_cost_m,sensor_cost_m,platform_cost_m,min_platforms_required,total_cost_m,norm_total_cost_m,'// &
                     'mission_max_reps,mission_max_reps_per_cost_m,'// &
                     'mission_reps_per_hour,mission_reps_per_hour_per_cost_m,'// &
                     'mission_one_rep_time_hr,cost_m_per_mission_one_rep_time_hr,'// &
                     'max_total_survey_area_nmi2,max_total_survey_area_nmi2_per_cost_m,'// &
                     'norm_max_reps,norm_max_reps_per_cost_m,norm_reps_per_hour,norm_reps_per_hour_per_cost_m,'// &
                     'norm_one_rep_time_hr,norm_cost_m_per_rep_time_hr,'// &
                     'norm_max_total_survey_area_nmi2,norm_max_total_survey_area_nmi2_per_cost_m,'// &
                     'average_absolute_score',absolute_performance_weight,',average_efficiency_score',efficiency_weight, &
                     ',average_total_score,rss_cost_one_rep_time,rss_norm_cost_norm_one_rep_time'
    write(fmt_str,'(a,i0,a)') '(e22.15,',size(results,dim=1)-1,'(",",e22.15))'
    do c_ii=1,size(results,dim=2)
        write(unit=fid, fmt=fmt_str) results(:,c_ii)
    end do
    close(fid)

end program main
