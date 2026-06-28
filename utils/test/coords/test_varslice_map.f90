program test_varslice_map
    ! Validate Stage A (SCRIP -> weight_map bridge) and the varslice migration to
    ! the unified map_field, against the legacy map_scrip_field applier.
    !
    !   1. Bridge: a hand-built conservative SCRIP map converted via
    !      map_scrip_to_weight_map and applied with map_field("mean") must match
    !      both the analytic expectation and the legacy map_scrip_field("mean").
    !   2. varslice: varslice_map_to_grid driven by the new map_class must produce
    !      the same target field and pick up its dimensions from map%G.

    use precision
    use mapping_scrip, only : map_scrip_class, map_scrip_field, map_scrip_to_weight_map
    use mapping,       only : map_class, map_field
    use varslice,      only : varslice_class, varslice_map_to_grid

    implicit none

    type(map_scrip_class) :: mps
    type(map_class)       :: map
    type(varslice_class)  :: vs_src, vs_tgt

    real(dp) :: var1d(2,2), v2_new(2,1), v2_leg(2,1), expected(2,1)
    real(wp) :: var1s(2,2)
    logical  :: m2(2,1)
    real(dp), parameter :: rtol = 1.0e-10_dp
    integer :: nfail

    nfail = 0

    ! --- hand-build a tiny conservative SCRIP map: 2x2 source -> 2x1 target ----
    ! dst 1 = 0.5*src1 + 0.5*src2 ; dst 2 = 0.5*src3 + 0.5*src4
    mps%src_name      = "src2x2"
    mps%dst_name      = "dst2x1"
    mps%src_grid_size = 4
    mps%dst_grid_size = 2
    mps%num_links     = 4
    mps%num_wgts      = 1
    allocate(mps%src_address(4), mps%dst_address(4), mps%remap_matrix(1,4))
    allocate(mps%dst_grid_dims(2))
    mps%src_address       = [1, 2, 3, 4]
    mps%dst_address       = [1, 1, 2, 2]      ! grouped by destination, ascending
    mps%remap_matrix(1,:) = [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp]
    mps%dst_grid_dims     = [2, 1]

    ! source field (column-major: src index 1..4) and analytic target
    var1d    = reshape([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp], [2,2])
    expected = reshape([15.0_dp, 35.0_dp], [2,1])

    ! --- Stage A: bridge + map_field("mean") vs legacy map_scrip_field("mean") --
    call map_scrip_to_weight_map(mps, map%wm)
    map%npts    = map%wm%n_dst
    map%is_grid = .true.
    map%G%nx    = 2
    map%G%ny    = 1

    v2_new = 0.0_dp
    call map_field(map, "v", var1d, v2_new, stat="mean", mask2=m2)

    v2_leg = 0.0_dp
    call map_scrip_field(mps, "v", var1d, v2_leg, method="mean")

    call check("bridge vs expected", maxval(abs(v2_new - expected)), rtol, nfail)
    call check("bridge vs legacy  ", maxval(abs(v2_new - v2_leg)),   rtol, nfail)
    if (.not. all(m2)) then
        write(*,*) "  FAIL mask2 not all true"
        nfail = nfail + 1
    end if

    ! --- varslice migration: varslice_map_to_grid on the new map_class ---------
    vs_src%par%name      = "v"
    vs_src%par%with_time = .false.
    allocate(vs_src%dim(2)); vs_src%dim = [2, 2]
    allocate(vs_src%x(2));   vs_src%x   = [0.0_wp, 1.0_wp]
    allocate(vs_src%y(2));   vs_src%y   = [0.0_wp, 1.0_wp]
    var1s = real(var1d, wp)
    allocate(vs_src%var(2,2,1,1)); vs_src%var(:,:,1,1) = var1s

    call varslice_map_to_grid(vs_tgt, vs_src, map, stat="mean")

    call check("varslice vs expected", &
               maxval(abs(real(vs_tgt%var(:,:,1,1),dp) - expected)), rtol, nfail)

    if (size(vs_tgt%var,1) /= 2 .or. size(vs_tgt%var,2) /= 1) then
        write(*,*) "  FAIL vs_tgt var dims wrong: ", shape(vs_tgt%var)
        nfail = nfail + 1
    end if

    if (nfail == 0) then
        write(*,*) "test_varslice_map: PASS"
    else
        write(*,*) "test_varslice_map: FAIL", nfail, "checks"
        stop 1
    end if

contains

    subroutine check(label, err, rtol, nfail)
        character(len=*), intent(in)    :: label
        real(dp),         intent(in)    :: err, rtol
        integer,          intent(inout) :: nfail
        if (err <= rtol) then
            write(*,"(a,a,es12.4)") "  PASS ", label//"  err=", err
        else
            write(*,"(a,a,es12.4)") "  FAIL ", label//"  err=", err
            nfail = nfail + 1
        end if
    end subroutine check

end program test_varslice_map
