program bench_apply
    ! Per-timestep APPLY throughput of map_field (NOT weight generation).
    !
    ! Build each map once, then time many repeated applies -- this is exactly
    ! what a coupled model (e.g. climber-x SMB) pays every step, where the map is
    ! cached and only map_field runs. Reports ms/apply and Mtarget-pts/s so a
    ! regression in the apply loop shows up directly.
    !
    ! The companion bench_alloc.f90 isolates *why* the apply loop is (or isn't)
    ! fast by A/B-timing per-target allocation vs a hoisted scratch on identical
    ! weight data.

    use coords
    implicit none

    integer, parameter :: NAPPLY = 30
    real(dp), parameter :: d2r = acos(-1.0_dp)/180.0_dp

    type(grid_class) :: g_glob, g_nh

    ! global 1-degree source
    call grid_init(g_glob, name="global-1deg", mtype="latlon", units="degrees", &
                   x0=0.0_dp, dx=1.0_dp, nx=360, y0=-90.0_dp, dy=1.0_dp, ny=181)

    ! NH polar-stereographic 16 km cap centered on the pole (+/-4000 km -> 500
    ! pts/side, ~250k targets). Representative of the SMB coupler->ice apply cost;
    ! apply time scales ~linearly with the target-point count.
    call grid_init(g_nh, name="NH-16KM", mtype="stereographic", units="kilometers", &
                   x0=-4000.0_dp, dx=16.0_dp, nx=500, &
                   y0=-4000.0_dp, dy=16.0_dp, ny=500, &
                   lambda=0.0_dp, phi=90.0_dp, alpha=71.0_dp)

    write(*,"(a)") ""
    write(*,"(a,i0,a,i0,a,i0)") "global src pts = ", g_glob%npts, &
        "    NH-16KM tgt pts = ", g_nh%npts, "    n_apply = ", NAPPLY
    write(*,"(a)") ""
    write(*,"(a)") "direction         method     n_links   t_init[s]  t_apply[ms]   Mtgt-pts/s"
    write(*,"(a)") "--------------------------------------------------------------------------"

    call bench_dir("global->NH-16KM", g_glob, g_nh, "con")
    call bench_dir("global->NH-16KM", g_glob, g_nh, "shepard")
    call bench_dir("NH-16KM->global", g_nh, g_glob, "con")
    call bench_dir("NH-16KM->global", g_nh, g_glob, "shepard")

    write(*,"(a)") ""
    write(*,"(a)") "done."

contains

    subroutine bench_dir(label, gs, gt, method)
        character(len=*), intent(in) :: label, method
        type(grid_class), intent(in) :: gs, gt

        type(map_class)       :: map
        real(dp), allocatable :: fs(:,:), ft(:,:)
        integer(8) :: c0, c1, cr
        integer    :: r
        real(dp)   :: t_init, t_apply

        allocate(fs(gs%G%nx, gs%G%ny), ft(gt%G%nx, gt%G%ny))
        if (gs%cs%is_cartesian) then
            fs = 1.0_dp + 1.0e-3_dp*(gs%x + gs%y)
        else
            fs = 2.0_dp + cos(gs%lat*d2r)*cos(gs%lon*d2r)
        end if

        ! ---- init (build weights once) ----
        call system_clock(c0, cr)
        if (trim(method) == "con") then
            call map_init_conservative(map, gs, gt)
        else
            call map_init(map, gs, gt, max_neighbors=10, load=.false.)
        end if
        call system_clock(c1)
        t_init = real(c1-c0,dp)/real(cr,dp)

        ! ---- apply (per-timestep hot path) ----
        call system_clock(c0, cr)
        do r = 1, NAPPLY
            if (trim(method) == "con") then
                call map_field(map, "f", fs, ft, stat="mean")
            else
                call map_field(map, "f", fs, ft, method="shepard")
            end if
        end do
        call system_clock(c1)
        t_apply = real(c1-c0,dp)/real(cr,dp)/real(NAPPLY,dp)

        write(*,"(a16,2x,a8,2x,i9,2x,es10.3,2x,f10.3,2x,f10.2)") &
            adjustl(label), adjustl(method), map%wm%n_links, t_init, &
            t_apply*1.0e3_dp, real(gt%npts,dp)/t_apply/1.0e6_dp
        flush(6)

        deallocate(fs, ft)
    end subroutine bench_dir

end program bench_apply
