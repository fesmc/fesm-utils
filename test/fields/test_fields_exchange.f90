program test_fields_exchange
    ! fesm_fields with coupler, as the kryos wiring layer will use them (the
    ! wire_t/exchange below are a sketch of kryos_wiring, design note §5):
    !   * con of an extensive field from a fine to a coarse grid conserves the
    !     area integral;
    !   * a rank-0 field bypasses the remap bit-exactly, validity forwarded;
    !   * MEAN -> MEAN with a consumer twice as slow as the producer: the
    !     import, finalized at the consumer's step, is the mean over the two
    !     producer intervals;
    !   * the remap workspace is not reallocated after the first call.
    ! Both grids cover [0,1000] km x [0,1000] km with aligned cell edges.

    use, intrinsic :: iso_c_binding, only : c_loc, c_associated, c_ptr
    use coords
    use coupler
    use fesm_log, only : fesm_error
    use fesm_fields

    implicit none

    type :: wire_t
        integer :: hs = 0, hd = 0
        character(len=16) :: method = ""
        real(dp), allocatable :: work(:,:,:)
    end type wire_t

    character(len=*), parameter :: fldr = "maps_fields_test"
    type(coupler_class) :: cpl
    type(grid_class)    :: gf, gc
    type(field_state_t) :: ice_ex, srf_im, srf_ex, ice_im
    type(wire_t), target :: w_dz, w_bsl, w_smb
    type(c_ptr) :: p_first
    integer  :: nfail, i, j, k, hdz, hbsl, hsmb
    real(dp) :: dz(20,20), dzc(10,10), bsl, bsl_in, t0, t1, smb(10,10), r(20,20,2), expect(20,20)
    real(dp) :: int_f, int_c

    nfail = 0
    call execute_command_line("rm -rf "//fldr)

    call grid_init(gf, name="fields-fine",   mtype="cartesian", units="kilometers", &
                   x0=25.0_dp, dx=50.0_dp,  nx=20, y0=25.0_dp, dy=50.0_dp,  ny=20)
    call grid_init(gc, name="fields-coarse", mtype="cartesian", units="kilometers", &
                   x0=50.0_dp, dx=100.0_dp, nx=10, y0=50.0_dp, dy=100.0_dp, ny=10)
    call coupler_init(cpl, map_fldr=fldr)
    call coupler_add_grid(cpl, "fields-fine",   gf)
    call coupler_add_grid(cpl, "fields-coarse", gc)

    ! ice (fine grid, dt = 2) <-> surface (coarse grid, dt = 1)
    call fields_init(ice_ex, "ice:exports")
    call fields_add(ice_ex, "dz",  "m", FIELD_INSTANT, FIELD_EXTENSIVE, hdz,  grid="fields-fine", dims=[20,20])
    call fields_add(ice_ex, "bsl", "m", FIELD_INSTANT, FIELD_INTENSIVE, hbsl)
    call fields_lock(ice_ex)

    call fields_init(srf_im, "srf:imports")
    call fields_add(srf_im, "dz",  "metres", FIELD_INSTANT, FIELD_EXTENSIVE, i, grid="fields-coarse", dims=[10,10])
    call fields_add(srf_im, "bsl", "m",      FIELD_INSTANT, FIELD_INTENSIVE, i)
    call fields_lock(srf_im)

    call fields_init(srf_ex, "srf:exports")
    call fields_add(srf_ex, "smb", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, hsmb, grid="fields-coarse", dims=[10,10])
    call fields_lock(srf_ex)

    call fields_init(ice_im, "ice:imports")
    call fields_add(ice_im, "smb", "mm w.e./yr", FIELD_MEAN, FIELD_EXTENSIVE, i, grid="fields-fine", dims=[20,20])
    call fields_lock(ice_im)

    call wire(w_dz,  ice_ex, "dz",  srf_im, "dz")
    call wire(w_bsl, ice_ex, "bsl", srf_im, "bsl")
    call wire(w_smb, srf_ex, "smb", ice_im, "smb")

    ! ----- con conserves the area integral of an extensive field ---------------
    do j = 1, 20
        do i = 1, 20
            dz(i,j) = 1.0_dp + sin(0.3_dp*real(i, dp))*cos(0.2_dp*real(j, dp)) + 0.01_dp*real(i*j, dp)
        end do
    end do
    bsl = -120.0_dp + 1.0_dp/3.0_dp
    call fields_put(ice_ex, hdz,  dz,  t1=2.0_dp)
    call fields_put(ice_ex, hbsl, bsl, t0=0.0_dp, t1=2.0_dp)

    call exchange(ice_ex, w_dz,  srf_im, cpl)
    call exchange(ice_ex, w_bsl, srf_im, cpl)

    call fields_get(srf_im, w_dz%hd, dzc)
    int_f = sum(dz)  * 50.0_dp**2
    int_c = sum(dzc) * 100.0_dp**2
    write(*,"(a,es12.4)") "  con: relative change of the area integral = ", abs(int_c - int_f)/abs(int_f)
    call check("con conserves the area integral (rel 1e-12)", abs(int_c - int_f) <= 1.0e-12_dp*abs(int_f))

    ! ----- rank 0 bypasses the remap --------------------------------------------
    call fields_get(srf_im, w_bsl%hd, bsl_in, t0=t0, t1=t1)
    call check("rank 0 bit-exact", bsl_in == bsl)
    call check("rank 0 validity forwarded", t0 == 0.0_dp .and. t1 == 2.0_dp)

    ! ----- MEAN -> MEAN, consumer slower than producer --------------------------
    ! The surface steps over [0,1] and [1,2]; after each step it finalizes its
    ! exports and the exchange extends the ice import window. The ice (dt = 2)
    ! finalizes its imports at t = 2.
    do k = 1, 2
        do j = 1, 10
            do i = 1, 10
                smb(i,j) = real(k, dp)*100.0_dp + real(i - j, dp)
            end do
        end do
        call fields_put(srf_ex, hsmb, smb, t0=real(k-1, dp), t1=real(k, dp))
        call fields_finalize_interval(srf_ex, real(k, dp))
        call exchange(srf_ex, w_smb, ice_im, cpl)
        if (k == 1) p_first = c_loc(w_smb%work(1,1,1))
        r(:,:,k) = w_smb%work(:,:,1)
    end do
    call check("workspace not reallocated after the first call", c_associated(p_first, c_loc(w_smb%work(1,1,1))))

    call fields_finalize_interval(ice_im, 2.0_dp)
    call fields_get(ice_im, w_smb%hd, expect, t0=t0, t1=t1)
    call check("slow consumer: mean over two producer intervals", &
               all(expect == (r(:,:,1)*1.0_dp + r(:,:,2)*1.0_dp)/2.0_dp))
    call check("slow consumer: validity [0,2]", t0 == 0.0_dp .and. t1 == 2.0_dp)

    call execute_command_line("rm -rf "//fldr)

    if (nfail > 0) then
        write(*,*) "test_fields_exchange: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_exchange"

contains

    subroutine wire(w, src, name_s, dst, name_d)
        ! Resolve a wiring row at init and check the link.
        type(wire_t),        intent(inout) :: w
        type(field_state_t), intent(in)    :: src, dst
        character(len=*),    intent(in)    :: name_s, name_d

        w%hs = fields_handle(src, name_s)
        w%hd = fields_handle(dst, name_d)
        call fields_check_link(src, w%hs, dst, w%hd)

    end subroutine wire

    subroutine exchange(src, w, dst, cpl)
        ! One wiring row: the published source value into the destination.
        type(field_state_t), intent(in)    :: src
        type(wire_t),        intent(inout) :: w
        type(field_state_t), intent(inout) :: dst
        type(coupler_class), intent(inout) :: cpl

        character(len=16) :: method

        associate (fs => src%f(w%hs))
            if (.not. fs%published) call fesm_error("exchange", "source not published", trim(fs%name))
            if (fs%rank == 0) then
                call fields_put(dst, w%hd, fs%v(1,1,1), t0=fs%t0, t1=fs%t1)
            else
                method = w%method
                if (len_trim(method) == 0) method = merge("con  ", "bilin", fs%cons == FIELD_EXTENSIVE)
                call remap(cpl, fs%v, trim(fs%grid), w%work, trim(dst%f(w%hd)%grid), method=trim(method))
                call fields_put(dst, w%hd, w%work, t0=fs%t0, t1=fs%t1)
            end if
        end associate

    end subroutine exchange

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

end program test_fields_exchange
