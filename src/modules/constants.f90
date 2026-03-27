module constants
use, non_intrinsic :: kinds, only: i32, sp, dp, c_bool
implicit none
private

#ifdef __COMPILE_FOR_DEBUG__
    logical(kind=c_bool), parameter :: debug = .true.
#else
    logical(kind=c_bool), parameter :: debug = .false.
#endif
    
    integer, parameter :: vec_len = 64_i32

    real(kind=sp), parameter :: eps_sp = sqrt(sqrt(spacing(0.0_sp)))
    real(kind=dp), parameter :: eps_dp = sqrt(sqrt(spacing(0.0_dp)))

    real(kind=sp), parameter :: pi_sp = acos(-1.0_sp)
    real(kind=dp), parameter :: pi_dp = acos(-1.0_dp)

    real(kind=sp), parameter :: deg2rad_sp = pi_sp/180.0_sp
    real(kind=dp), parameter :: deg2rad_dp = pi_dp/180.0_dp
    real(kind=sp), parameter :: rad2deg_sp = 180.0_sp/pi_sp
    real(kind=dp), parameter :: rad2deg_dp = 180.0_dp/pi_dp

    real(kind=sp), parameter :: twopi_sp = 2.0_sp*acos(-1.0_sp)
    real(kind=dp), parameter :: twopi_dp = 2.0_dp*acos(-1.0_dp)

    real(kind=sp), parameter :: nmi2ft_sp = 1852.0_sp*100.0_sp/2.54_sp/12.0_sp
    real(kind=dp), parameter :: nmi2ft_dp = 1852.0_dp*100.0_dp/2.54_dp/12.0_dp

    real(kind=sp), parameter :: ft2nmi_sp = 2.54_sp*12.0_sp/1852.0_sp/100.0_sp
    real(kind=dp), parameter :: ft2nmi_dp = 2.54_dp*12.0_dp/1852.0_dp/100.0_dp

    real(kind=sp), parameter :: km2nmi_sp = 1.0_sp/1.852_sp
    real(kind=dp), parameter :: km2nmi_dp = 1.0_dp/1.852_dp

    public :: debug, pi_sp, pi_dp, vec_len, twopi_sp, twopi_dp, deg2rad_sp, deg2rad_dp, rad2deg_sp, rad2deg_dp, eps_sp, eps_dp, &
              nmi2ft_sp, nmi2ft_dp, ft2nmi_sp, ft2nmi_dp, km2nmi_sp, km2nmi_dp

end module constants
