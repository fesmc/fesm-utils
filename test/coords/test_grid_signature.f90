program test_grid_signature
    ! Map cache identity (docs/dev/coupler-global-design.md, 4):
    !   * a grid and its grid_<name>.txt round trip have matching signatures
    !     (and equal hash tags); a grid with another dy does not match,
    !   * two grids with one name and different dy in one folder:
    !     on_mismatch="error" stops (child process, exit code checked),
    !     on_mismatch="rebuild" regenerates the file,
    !   * a legacy cache file without signatures loads (unverified),
    !   * a cache file of the wrong size stops (child process),
    !   * hash_names=.true. gives two files that coexist.
    ! Run without arguments. The child modes ("mismatch", "size") are invoked
    ! by the program itself.

    use coords

    implicit none

    character(len=*), parameter :: fldr = "maps_sig_test"
    type(grid_class) :: ga, gb, gb2, gr
    type(map_class)  :: map
    character(len=512) :: self, arg, why
    character(len=8)   :: h_a, h_r
    integer :: fails, nfiles, rc
    real(dp), allocatable :: f1(:,:), f2(:,:)

    call get_command_argument(0, self)
    arg = ""
    if (command_argument_count() > 0) call get_command_argument(1, arg)

    call make_grids(ga, gb, gb2)

    select case (trim(arg))
        case ("mismatch")
            ! the file for (sig-a -> sig-b) was built with gb; gb2 has another dy
            call map_init(map, ga, gb2, method="con", fldr=fldr)
            stop 0          ! not reached: map_init must stop
        case ("size")
            call map_init(map, ga, gb, method="nn", fldr=fldr)
            stop 0          ! not reached
    end select

    fails = 0
    call execute_command_line("rm -rf "//fldr//" ; mkdir -p "//fldr)

    ! --- 1. signatures: round trip matches, another dy does not ------------------
    call grid_cdo_write_desc_short(ga, fldr)
    call grid_cdo_read_desc(gr, trim(ga%name), fldr)
    if (.not. grid_signature_match(grid_signature_str(ga), grid_signature_str(gr), why)) then
        write(*,*) "FAIL: round trip signature differs: ", trim(why); fails = fails + 1
    end if
    h_a = grid_hash8(ga); h_r = grid_hash8(gr)
    if (h_a /= h_r) then
        write(*,*) "FAIL: round trip hash differs: ", h_a, " ", h_r; fails = fails + 1
    end if
    if (grid_signature_match(grid_signature_str(gb), grid_signature_str(gb2), why)) then
        write(*,*) "FAIL: grids with different dy have matching signatures"; fails = fails + 1
    else
        write(*,*) "signature: round trip matches (hash ", h_a, "), other dy differs: ", trim(why)
    end if
    if (grid_hash8(gb) == grid_hash8(gb2)) then
        write(*,*) "FAIL: grids with different dy have one hash"; fails = fails + 1
    end if

    ! --- 2. same name, other grid: error (child) / rebuild ----------------------------
    call map_init(map, ga, gb, method="con", fldr=fldr)                 ! writes the signed file
    call execute_command_line(trim(self)//" mismatch > "//fldr//"/child_mismatch.log 2>&1", exitstat=rc)
    if (rc == 0 .or. .not. log_has(fldr//"/child_mismatch.log", "built for another target grid")) then
        write(*,*) "FAIL: on_mismatch='error' did not stop on the signature"; fails = fails + 1
    else
        write(*,*) "on_mismatch='error' stops (exit ", rc, ")"
    end if
    call map_init(map, ga, gb2, method="con", fldr=fldr, on_mismatch="rebuild")
    if (map%wm%n_dst /= size(gb2%x)) then
        write(*,*) "FAIL: rebuild did not regenerate for the new grid"; fails = fails + 1
    end if
    ! the rebuilt file now matches gb2 and loads without a rebuild
    call map_init(map, ga, gb2, method="con", fldr=fldr)
    allocate(f1(size(ga%x,1), size(ga%x,2)), f2(size(gb2%x,1), size(gb2%x,2)))
    f1 = 1.0_dp
    call map_field(map, "one", f1, f2)
    if (abs(maxval(f2) - 1.0_dp) > 1.0e-12_dp) then
        write(*,*) "FAIL: rebuilt map does not preserve a constant"; fails = fails + 1
    else
        write(*,*) "on_mismatch='rebuild' regenerates"
    end if

    ! --- 3. legacy file (no signatures) loads ------------------------------------------
    call map_init(map, ga, gb, method="bilin", fldr=fldr, precompute_weights=.false.)
    call weight_map_write(map%wm, fldr//"/scrip-shepard_"//trim(ga%name)//"_"//trim(gb%name)//".nc")
    call map_init(map, ga, gb, method="shepard", fldr=fldr)
    if (map%wm%n_dst /= size(gb%x)) then
        write(*,*) "FAIL: legacy file did not load"; fails = fails + 1
    else
        write(*,*) "legacy file loads (unverified)"
    end if

    ! --- 4. wrong size stops (child) ---------------------------------------------------
    ! a legacy nn file for (ga -> ga) put under the (ga -> gb) name
    call map_init(map, ga, ga, method="nn", fldr=fldr, precompute_weights=.false.)
    call weight_map_write(map%wm, fldr//"/scrip-nn_"//trim(ga%name)//"_"//trim(gb%name)//".nc")
    call execute_command_line(trim(self)//" size > "//fldr//"/child_size.log 2>&1", exitstat=rc)
    if (rc == 0 .or. .not. log_has(fldr//"/child_size.log", "does not match the grid sizes")) then
        write(*,*) "FAIL: wrong-size map file did not stop on the size check"; fails = fails + 1
    else
        write(*,*) "wrong-size map file stops (exit ", rc, ")"
    end if

    ! --- 5. hash_names: two grids with one name coexist ----------------------------
    call execute_command_line("rm -f "//fldr//"/scrip-con_*")
    call map_init(map, ga, gb,  method="con", fldr=fldr, hash_names=.true.)
    call map_init(map, ga, gb2, method="con", fldr=fldr, hash_names=.true.)
    call map_init(map, ga, gb,  method="con", fldr=fldr, hash_names=.true.)   ! loads, no error
    call execute_command_line("ls "//fldr//"/scrip-con_*.nc | wc -l > "//fldr//"/count.txt")
    open(10, file=fldr//"/count.txt"); read(10,*) nfiles; close(10)
    if (nfiles /= 2 .or. map%wm%n_dst /= size(gb%x)) then
        write(*,*) "FAIL: hash_names: expected 2 coexisting files, found", nfiles; fails = fails + 1
    else
        write(*,*) "hash_names: 2 files coexist"
    end if

    call execute_command_line("rm -rf "//fldr)
    if (fails > 0) then
        write(*,*) "FAIL: test_grid_signature (", fails, " failure(s))"
        stop 1
    end if
    write(*,*) "PASS: test_grid_signature"

contains

    logical function log_has(file, text)
        ! Does the text file contain `text`?
        character(len=*), intent(in) :: file, text
        integer :: u, ios
        character(len=1024) :: line
        log_has = .false.
        open(newunit=u, file=file, status="old", action="read", iostat=ios)
        if (ios /= 0) return
        do
            read(u, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (index(line, text) > 0) log_has = .true.
        end do
        close(u)
    end function log_has

    subroutine make_grids(ga, gb, gb2)
        type(grid_class), intent(out) :: ga, gb, gb2
        call grid_init(ga,  name="sig-a", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                       lon180=.false., x0=-400.0_dp, dx=40.0_dp, nx=21, y0=-400.0_dp, dy=40.0_dp, ny=21, &
                       lambda=0.0_dp, phi=-71.0_dp)
        call grid_init(gb,  name="sig-b", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                       lon180=.false., x0=-400.0_dp, dx=80.0_dp, nx=11, y0=-400.0_dp, dy=80.0_dp, ny=11, &
                       lambda=0.0_dp, phi=-71.0_dp)
        ! same name and size as gb, other dy
        call grid_init(gb2, name="sig-b", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                       lon180=.false., x0=-400.0_dp, dx=80.0_dp, nx=11, y0=-400.0_dp, dy=70.0_dp, ny=11, &
                       lambda=0.0_dp, phi=-71.0_dp)
    end subroutine make_grids

end program test_grid_signature
