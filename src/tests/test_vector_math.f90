program test_vector_math
use, non_intrinsic :: kinds, only: stdout, sp, dp
use, non_intrinsic :: vector_math, only: vmag, vunit
implicit none

    real(kind=sp) :: v2sp(2), v3sp(3), vsp(4), unitsp(4)
    real(kind=dp) :: v2dp(2), v3dp(3), vdp(4), unitdp(4)

    v2sp = [3.0_sp, 4.0_sp]
    call vunit(v2sp, unitsp(1:2))
    write(unit=stdout, fmt='(a,2f4.1,a,f0.1,a,2f4.1)') 'v2sp: ',v2sp,', magnitude: ',vmag(v2sp),', unit: ',unitsp(1:2)

    v3sp = [3.0_sp, 0.0_sp, 4.0_sp]
    call vunit(v3sp, unitsp(1:3))
    write(unit=stdout, fmt='(a,3f4.1,a,f0.1,a,3f4.1)') 'v3sp: ',v3sp,', magnitude: ',vmag(v3sp),', unit: ',unitsp(1:3)

    vsp = [0.0_sp, 3.0_sp, 0.0_sp, 4.0_sp]
    call vunit(vsp, unitsp)
    write(unit=stdout, fmt='(a,4f4.1,a,f0.1,a,4f4.1)') 'vsp: ',vsp,', magnitude: ',vmag(vsp),', unit: ',unitsp
    write(unit=stdout, fmt='(3(a,f0.1))') 'vmag(v2sp): ',vmag(v2sp),', vmag(v3sp): ',vmag(v3sp),', vmag(vsp): ',vmag(vsp)

    v2dp = [3.0_dp, 4.0_dp]
    call vunit(v2dp, unitdp(1:2))
    write(unit=stdout, fmt='(a,2f4.1,a,f0.1,a,2f4.1)') 'v2dp: ',v2dp,', magnitude: ',vmag(v2dp),', unit: ',unitdp(1:2)

    v3dp = [3.0_dp, 0.0_dp, 4.0_dp]
    call vunit(v3dp, unitdp(1:3))
    write(unit=stdout, fmt='(a,3f4.1,a,f0.1,a,3f4.1)') 'v3dp: ',v3dp,', magnitude: ',vmag(v3dp),', unit: ',unitdp(1:3)

    vdp = [0.0_dp, 3.0_dp, 0.0_dp, 4.0_dp]
    call vunit(vdp, unitdp)
    write(unit=stdout, fmt='(a,4f4.1,a,f0.1,a,4f4.1)') 'vdp: ',vdp,', magnitude: ',vmag(vdp),', unit: ',unitdp
    write(unit=stdout, fmt='(3(a,f0.1))') 'vmag(v2dp): ',vmag(v2dp),', vmag(v3dp): ',vmag(v3dp),', vmag(vdp): ',vmag(vdp)

end program test_vector_math
