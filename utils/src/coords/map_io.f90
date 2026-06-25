module map_io
    ! Read/write a weight_map as a SCRIP-superset NetCDF file.
    !
    ! The file is always a valid SCRIP remapping file (src_address, dst_address,
    ! remap_matrix), so it interoperates with cdo and other tools. A MAP_DISTANCE
    ! map additionally writes coords-only extension variables (mp_dist,
    ! mp_quadrant, mp_xn, mp_yn) so it can be reloaded losslessly and keep
    ! change-method-for-free; the SCRIP core then holds shepard-baked weights so
    ! the file is still usable as a standard map. With scrip_only=.true. only the
    ! standard SCRIP core is written and a reload yields a frozen MAP_WEIGHT.

    use precision,   only: dp
    use weight_map,  only: weight_map_t, weight_map_alloc, weight_map_index, &
                           MAP_DISTANCE, MAP_WEIGHT
    use ncio

    implicit none
    private

    public :: weight_map_write, weight_map_read

contains

    subroutine weight_map_write(wm, filename, scrip_only)
        type(weight_map_t), intent(in) :: wm
        character(len=*),   intent(in) :: filename
        logical, optional,  intent(in) :: scrip_only

        logical  :: sonly
        integer  :: i, nl, store_kind, meta(3)
        real(dp), allocatable :: wbake(:)

        sonly = .false.
        if (present(scrip_only)) sonly = scrip_only
        nl = wm%n_links

        ! SCRIP-core weights
        allocate(wbake(nl))
        if (wm%kind == MAP_WEIGHT) then
            wbake = wm%w
        else
            call bake_shepard(wm, wbake)
        end if

        ! kind as stored: a scrip_only distance map reloads as a frozen weight map
        store_kind = wm%kind
        if (wm%kind == MAP_DISTANCE .and. sonly) store_kind = MAP_WEIGHT

        call nc_create(filename)
        call nc_write_dim(filename, "num_links", x=[(dble(i), i=1,nl)])
        call nc_write_dim(filename, "num_wgts",  x=[1.0_dp])
        call nc_write_dim(filename, "params",    x=[(dble(i), i=1,3)])

        meta = [wm%n_src, wm%n_dst, store_kind]
        call nc_write(filename, "meta", meta, dim1="params")
        call nc_write(filename, "src_address", wm%src, dim1="num_links")
        call nc_write(filename, "dst_address", wm%dst, dim1="num_links")
        call nc_write(filename, "remap_matrix", reshape(wbake, [1, nl]), &
                      dim1="num_wgts", dim2="num_links")

        if (wm%kind == MAP_DISTANCE .and. .not. sonly) then
            call nc_write(filename, "mp_dist",     wm%dist,     dim1="num_links")
            call nc_write(filename, "mp_quadrant", wm%quadrant, dim1="num_links")
            call nc_write(filename, "mp_xn",       wm%xn,       dim1="num_links")
            call nc_write(filename, "mp_yn",       wm%yn,       dim1="num_links")
        end if

        deallocate(wbake)
    end subroutine weight_map_write

    subroutine weight_map_read(wm, filename)
        type(weight_map_t), intent(inout) :: wm
        character(len=*),   intent(in)    :: filename
        integer :: nl, kind, n_src, n_dst, meta(3)
        character(len=56), allocatable :: dim_names(:)
        integer,           allocatable :: dims(:)
        real(dp), allocatable :: rmat(:,:)

        call nc_dims(filename, "src_address", dim_names, dims)
        nl = dims(1)
        call nc_read(filename, "meta", meta)
        n_src = meta(1); n_dst = meta(2); kind = meta(3)

        call weight_map_alloc(wm, kind, n_src, n_dst, nl)
        call nc_read(filename, "src_address", wm%src)
        call nc_read(filename, "dst_address", wm%dst)

        if (kind == MAP_WEIGHT) then
            allocate(rmat(1, nl))
            call nc_read(filename, "remap_matrix", rmat)
            wm%w = rmat(1, :)
            deallocate(rmat)
        else
            call nc_read(filename, "mp_dist",     wm%dist)
            call nc_read(filename, "mp_quadrant", wm%quadrant)
            call nc_read(filename, "mp_xn",       wm%xn)
            call nc_read(filename, "mp_yn",       wm%yn)
        end if

        call weight_map_index(wm)
    end subroutine weight_map_read

    subroutine bake_shepard(wm, wbake)
        ! Per-target normalized inverse-distance-squared weights (sum to 1 per
        ! target), so the SCRIP core of a distance map is a usable shepard map.
        type(weight_map_t), intent(in)  :: wm
        real(dp),           intent(out) :: wbake(:)
        integer  :: k, j, j1, j2
        real(dp) :: wsum
        wbake = 0.0_dp
        do k = 1, wm%n_dst
            j1 = wm%dst_off(k); j2 = wm%dst_off(k+1) - 1
            if (j2 < j1) cycle
            wsum = 0.0_dp
            do j = j1, j2
                wbake(j) = 1.0_dp / max(wm%dist(j), 1.0e-12_dp)**2
                wsum = wsum + wbake(j)
            end do
            if (wsum > 0.0_dp) wbake(j1:j2) = wbake(j1:j2) / wsum
        end do
    end subroutine bake_shepard

end module map_io
