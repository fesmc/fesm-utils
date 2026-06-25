module conservative
    ! Analytic conservative remapping via polygon clipping (improvement-doc §2,
    ! stage B: planar Sutherland-Hodgman). For each target cell, candidate source
    ! cells are found with a k-d-tree radius query, then each source cell polygon
    ! is clipped against the target cell and the overlap area becomes the link
    ! weight. The result is a MAP_WEIGHT applied as an area-weighted mean.
    !
    ! Two planar modes are supported:
    !   - same-system: source and target share a Cartesian/projected plane;
    !     cells are clipped directly.
    !   - lat-lon source -> projected target: source cell corners (lon/lat) are
    !     projected into the target projection plane, then clipped there.
    ! Lat-lon / Gaussian *targets* need great-circle clipping (stage C, later).

    use precision,   only: dp
    use coordinates, only: grid_class, compare_coord
    use oblimap_projection_module, only: oblimap_projection
    use kdtree
    use weight_map
    use polygons
    use mapping,     only: map_class

    implicit none
    private

    public :: map_init_conservative

contains

    subroutine map_init_conservative(map, grid1, grid2)
        type(map_class),  intent(inout) :: map
        type(grid_class), intent(in)    :: grid1, grid2

        type(kdtree_t) :: tree
        logical  :: same_sys, do_project
        integer  :: n1, n2, nx1, ny1, nx2, ny2
        integer  :: i, j, ic, jc, c, cs, cc, k, nf, nl, nlmax, no
        real(dp) :: rad, diag2, area, maxsrcdiag, d, xyc, px, py
        real(dp) :: lx(4), ly(4), tcx(4), tcy(4), qpt(2)
        real(dp), allocatable :: emb(:,:), ox(:), oy(:)
        real(dp), allocatable :: scornx(:,:), scorny(:,:)   ! source corners in target plane
        integer,  allocatable :: cidx(:)
        real(dp), allocatable :: cd2(:)
        integer,  allocatable :: tsrc(:), tdst(:)
        real(dp), allocatable :: tw(:)

        same_sys   = compare_coord(grid1, grid2) .and. grid2%cs%is_cartesian
        do_project = (.not. same_sys) .and. grid2%cs%is_projection .and. (.not. grid1%cs%is_cartesian)

        if (.not. (same_sys .or. do_project)) then
            write(*,*) "conservative:: map_init_conservative: stage-B planar clip supports"
            write(*,*) "  (a) same Cartesian/projected system, or"
            write(*,*) "  (b) lat-lon source -> projected target."
            write(*,*) "  Lat-lon/Gaussian targets need the spherical clip (stage C)."
            stop
        end if

        nx1 = grid1%G%nx; ny1 = grid1%G%ny; n1 = nx1*ny1
        nx2 = grid2%G%nx; ny2 = grid2%G%ny; n2 = nx2*ny2
        xyc = grid2%cs%xy_conv

        map%name1 = grid1%name
        map%name2 = grid2%name
        map%cs    = grid2%cs
        map%is_grid = .true.
        map%G     = grid2%G
        map%npts  = n2
        map%nmax  = 0
        map%is_same_map = same_sys
        if (allocated(map%x)) deallocate(map%x, map%y, map%lon, map%lat)
        allocate(map%x(n2), map%y(n2), map%lon(n2), map%lat(n2))
        map%x   = reshape(grid2%x,   [n2])
        map%y   = reshape(grid2%y,   [n2])
        map%lon = reshape(grid2%lon, [n2])
        map%lat = reshape(grid2%lat, [n2])

        ! Express every source cell (corners + center) in the TARGET plane.
        allocate(scornx(4, n1), scorny(4, n1), emb(2, n1))
        maxsrcdiag = 0.0_dp
        do jc = 1, ny1
            do ic = 1, nx1
                cs = (jc-1)*nx1 + ic
                call cell_corners_grid(grid1, ic, jc, lx, ly)
                if (same_sys) then
                    scornx(:,cs) = lx; scorny(:,cs) = ly
                    emb(1,cs) = grid1%x(ic,jc); emb(2,cs) = grid1%y(ic,jc)
                else
                    ! lx,ly are lon,lat -> project into the target projection plane
                    do k = 1, 4
                        call oblimap_projection(lx(k), ly(k), px, py, grid2%cs%proj)
                        scornx(k,cs) = px/xyc; scorny(k,cs) = py/xyc
                    end do
                    call oblimap_projection(grid1%lon(ic,jc), grid1%lat(ic,jc), px, py, grid2%cs%proj)
                    emb(1,cs) = px/xyc; emb(2,cs) = py/xyc
                end if
                d = sqrt((maxval(scornx(:,cs))-minval(scornx(:,cs)))**2 &
                       + (maxval(scorny(:,cs))-minval(scorny(:,cs)))**2)
                maxsrcdiag = max(maxsrcdiag, d)
            end do
        end do
        call kdtree_build(tree, emb)

        ! candidate radius: half target diagonal + largest source-cell diagonal
        diag2 = sqrt(grid2%G%dx**2 + grid2%G%dy**2)
        rad   = 0.5_dp*diag2 + maxsrcdiag

        allocate(cidx(n1), cd2(n1))
        nlmax = max(n2*16, n1)
        allocate(tsrc(nlmax), tdst(nlmax), tw(nlmax))

        nl = 0
        do j = 1, ny2
            do i = 1, nx2
                c = (j-1)*nx2 + i                 ! target cell (column-major)
                call cell_corners_grid(grid2, i, j, tcx, tcy)
                qpt = [grid2%x(i,j), grid2%y(i,j)]
                call kdtree_radius(tree, qpt, rad, cidx, cd2, nf)
                do cc = 1, nf
                    cs = cidx(cc)
                    call polygon_clip(scornx(:,cs), scorny(:,cs), tcx, tcy, ox, oy, no)
                    if (no >= 3) then
                        area = polygon_area(ox(1:no), oy(1:no))
                        if (area > 0.0_dp) then
                            if (nl >= nlmax) call grow_links(tsrc, tdst, tw, nlmax)
                            nl = nl + 1
                            tsrc(nl) = cs
                            tdst(nl) = c
                            tw(nl)   = area
                        end if
                    end if
                end do
            end do
        end do

        call weight_map_alloc(map%wm, MAP_WEIGHT, n1, n2, nl)
        map%wm%src(1:nl) = tsrc(1:nl)
        map%wm%dst(1:nl) = tdst(1:nl)
        map%wm%w(1:nl)   = tw(1:nl)
        call weight_map_index(map%wm)

        call kdtree_free(tree)
        deallocate(emb, scornx, scorny, cidx, cd2, tsrc, tdst, tw)
        if (allocated(ox)) deallocate(ox, oy)
    end subroutine map_init_conservative

    subroutine grow_links(tsrc, tdst, tw, cap)
        ! Double the temporary link buffers (safety for dense overlaps).
        integer,  allocatable, intent(inout) :: tsrc(:), tdst(:)
        real(dp), allocatable, intent(inout) :: tw(:)
        integer,               intent(inout) :: cap
        integer,  allocatable :: itmp(:)
        real(dp), allocatable :: rtmp(:)
        integer :: newcap
        newcap = cap*2
        allocate(itmp(newcap)); itmp(1:cap) = tsrc; call move_alloc(itmp, tsrc)
        allocate(itmp(newcap)); itmp(1:cap) = tdst; call move_alloc(itmp, tdst)
        allocate(rtmp(newcap)); rtmp(1:cap) = tw;   call move_alloc(rtmp, tw)
        cap = newcap
    end subroutine grow_links

    subroutine cell_corners_grid(grid, i, j, cx, cy)
        ! Four corners (CCW) of cell (i,j) in axis units, from adjacent-axis
        ! midpoints (so uniform and non-uniform/Gaussian grids share one path).
        type(grid_class), intent(in)  :: grid
        integer,          intent(in)  :: i, j
        real(dp),         intent(out) :: cx(4), cy(4)
        real(dp) :: x1, x2, y1, y2
        integer  :: nx, ny
        nx = grid%G%nx; ny = grid%G%ny

        if (i == 1) then
            x1 = grid%G%x(i) - 0.5_dp*(grid%G%x(i+1)-grid%G%x(i))
            x2 = 0.5_dp*(grid%G%x(i)+grid%G%x(i+1))
        else if (i == nx) then
            x1 = 0.5_dp*(grid%G%x(i)+grid%G%x(i-1))
            x2 = grid%G%x(i) + 0.5_dp*(grid%G%x(i)-grid%G%x(i-1))
        else
            x1 = 0.5_dp*(grid%G%x(i)+grid%G%x(i-1))
            x2 = 0.5_dp*(grid%G%x(i)+grid%G%x(i+1))
        end if

        if (j == 1) then
            y1 = grid%G%y(j) - 0.5_dp*(grid%G%y(j+1)-grid%G%y(j))
            y2 = 0.5_dp*(grid%G%y(j)+grid%G%y(j+1))
        else if (j == ny) then
            y1 = 0.5_dp*(grid%G%y(j)+grid%G%y(j-1))
            y2 = grid%G%y(j) + 0.5_dp*(grid%G%y(j)-grid%G%y(j-1))
        else
            y1 = 0.5_dp*(grid%G%y(j)+grid%G%y(j-1))
            y2 = 0.5_dp*(grid%G%y(j)+grid%G%y(j+1))
        end if

        cx = [x1, x2, x2, x1]
        cy = [y1, y1, y2, y2]
    end subroutine cell_corners_grid

end module conservative
