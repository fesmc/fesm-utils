program test_ncio_interp
    ! Validate the migrated nc_read_interp family (now in module ncio_interp,
    ! applying via the unified map_field). The source field is supplied through the
    ! var_in argument so no NetCDF file is needed; the file-read branch is a thin
    ! ncio wrapper around the same 2D apply exercised here.
    !
    !   - copy path  : source dims == target dims -> pass-through.
    !   - remap path : source dims /= target dims -> remap with a map_class, for
    !                  dp / sp / integer / logical fields.

    use precision
    use mapping_scrip, only : map_scrip_class, map_scrip_to_weight_map
    use mapping,       only : map_class
    use ncio_interp,   only : nc_read_interp

    implicit none

    type(map_scrip_class) :: mps
    type(map_class)       :: map

    real(dp) :: src_dp(2,2), tgt_dp(2,1), copy_dp(2,2), exp_dp(2,1)
    real(sp) :: src_sp(2,2), tgt_sp(2,1)
    integer  :: src_in(2,2), tgt_in(2,1)
    logical  :: src_lg(2,2), tgt_lg(2,1)
    real(dp), parameter :: rtol = 1.0e-10_dp
    integer :: nfail

    nfail = 0

    ! --- hand-build a conservative map: 2x2 source -> 2x1 target ---------------
    ! dst 1 = 0.5*src1 + 0.5*src2 ; dst 2 = 0.5*src3 + 0.5*src4
    mps%src_grid_size = 4
    mps%dst_grid_size = 2
    mps%num_links     = 4
    mps%num_wgts      = 1
    allocate(mps%src_address(4), mps%dst_address(4), mps%remap_matrix(1,4), mps%dst_grid_dims(2))
    mps%src_address       = [1, 2, 3, 4]
    mps%dst_address       = [1, 1, 2, 2]
    mps%remap_matrix(1,:) = [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp]
    mps%dst_grid_dims     = [2, 1]

    call map_scrip_to_weight_map(mps, map%wm)
    map%npts    = map%wm%n_dst
    map%is_grid = .true.
    map%G%nx    = 2
    map%G%ny    = 1

    src_dp = reshape([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp], [2,2])
    exp_dp = reshape([15.0_dp, 35.0_dp], [2,1])

    ! --- copy path (dims match, no map needed) --------------------------------
    copy_dp = -1.0_dp
    call nc_read_interp("none", "v", copy_dp, src_dp)
    call check("dp copy path", maxval(abs(copy_dp - src_dp)), rtol, nfail)

    ! --- remap path, dp -------------------------------------------------------
    tgt_dp = 0.0_dp
    call nc_read_interp("none", "v", tgt_dp, src_dp, map=map, stat="mean")
    call check("dp remap     ", maxval(abs(tgt_dp - exp_dp)), rtol, nfail)

    ! --- remap path, sp -------------------------------------------------------
    src_sp = real(src_dp, sp)
    tgt_sp = 0.0_sp
    call nc_read_interp("none", "v", tgt_sp, src_sp, map=map, stat="mean")
    call check("sp remap     ", maxval(abs(real(tgt_sp,dp) - exp_dp)), rtol, nfail)

    ! --- remap path, integer --------------------------------------------------
    src_in = nint(src_dp)
    tgt_in = 0
    call nc_read_interp("none", "v", tgt_in, src_in, map=map, stat="mean")
    call check("int remap    ", maxval(abs(real(tgt_in,dp) - exp_dp)), rtol, nfail)

    ! --- remap path, logical (count + nn fill) --------------------------------
    ! src region for dst1 all true, for dst2 all false -> expect [T,F]
    src_lg = reshape([.true., .true., .false., .false.], [2,2])
    call nc_read_interp("none", "v", tgt_lg, src_lg, map=map)
    if (tgt_lg(1,1) .and. (.not. tgt_lg(2,1))) then
        write(*,"(a)") "  PASS logical remap  -> [T,F]"
    else
        write(*,"(a,2l2)") "  FAIL logical remap  got ", tgt_lg(1,1), tgt_lg(2,1)
        nfail = nfail + 1
    end if

    if (nfail == 0) then
        write(*,*) "test_ncio_interp: PASS"
    else
        write(*,*) "test_ncio_interp: FAIL", nfail, "checks"
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

end program test_ncio_interp
