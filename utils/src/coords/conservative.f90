module conservative
    ! Analytic conservative remapping via polygon clipping (improvement-doc §2,
    ! stage B: planar Sutherland-Hodgman). For each target cell, candidate source
    ! cells are found with a k-d-tree radius query, then each source cell polygon
    ! is clipped against the target cell and the overlap area becomes the link
    ! weight. The result is a MAP_WEIGHT applied as an area-weighted mean.
    !
    ! This stage handles target grids whose cells live in a shared planar
    ! (Cartesian / same-projection) system as the source. Cross-system and
    ! spherical (great-circle) clipping are later stages.

    use precision,   only: dp
    use coordinates, only: grid_class, compare_coord
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
        integer  :: n1, n2, nx1, ny1, nx2, ny2
        integer  :: i, j, ic, jc, c, cc, nf, nl, nlmax, no
        real(dp) :: rad, diag1, diag2, area
        real(dp) :: scx(4), scy(4), tcx(4), tcy(4), qpt(2)
        real(dp), allocatable :: emb(:,:), ox(:), oy(:)
        integer,  allocatable :: cidx(:)
        real(dp), allocatable :: cd2(:)
        integer,  allocatable :: tsrc(:), tdst(:)
        real(dp), allocatable :: tw(:)

        if (.not. (compare_coord(grid1, grid2) .and. grid2%cs%is_cartesian)) then
            write(*,*) "conservative:: map_init_conservative: stage-B planar clip requires"
            write(*,*) "  source and target to share a Cartesian/projected system."
            write(*,*) "  (cross-system and spherical clipping are later stages.)"
            stop
        end if

        nx1 = grid1%G%nx; ny1 = grid1%G%ny; n1 = nx1*ny1
        nx2 = grid2%G%nx; ny2 = grid2%G%ny; n2 = nx2*ny2

        map%name1 = grid1%name
        map%name2 = grid2%name
        map%cs    = grid2%cs
        map%is_grid = .true.
        map%G     = grid2%G
        map%npts  = n2
        map%nmax  = 0
        map%is_same_map = .true.
        if (allocated(map%x)) deallocate(map%x, map%y, map%lon, map%lat)
        allocate(map%x(n2), map%y(n2), map%lon(n2), map%lat(n2))
        map%x   = reshape(grid2%x,   [n2])
        map%y   = reshape(grid2%y,   [n2])
        map%lon = reshape(grid2%lon, [n2])
        map%lat = reshape(grid2%lat, [n2])

        ! k-d tree over source cell centers (planar x,y in axis units)
        allocate(emb(2, n1))
        emb(1,:) = reshape(grid1%x, [n1])
        emb(2,:) = reshape(grid1%y, [n1])
        call kdtree_build(tree, emb)

        ! candidate search radius: target diag + source diag (generous)
        diag1 = sqrt(grid1%G%dx**2 + grid1%G%dy**2)
        diag2 = sqrt(grid2%G%dx**2 + grid2%G%dy**2)
        rad   = 0.5_dp*diag1 + 0.5_dp*diag2 + 0.25_dp*max(diag1, diag2)

        allocate(cidx(n1), cd2(n1))
        nlmax = 0
        ! generous link estimate: each target overlaps a handful of source cells
        nlmax = n2 * max(16, int(2.0_dp*(diag2/max(diag1,1.0e-12_dp) + 1.0_dp))**2)
        nlmax = max(nlmax, n2*8)
        allocate(tsrc(nlmax), tdst(nlmax), tw(nlmax))

        nl = 0
        do j = 1, ny2
            do i = 1, nx2
                c = (j-1)*nx2 + i                 ! target cell (column-major)
                call cell_corners_grid(grid2, i, j, tcx, tcy)
                qpt = [grid2%x(i,j), grid2%y(i,j)]
                call kdtree_radius(tree, qpt, rad, cidx, cd2, nf)
                do cc = 1, nf
                    ic = mod(cidx(cc)-1, nx1) + 1
                    jc = (cidx(cc)-1)/nx1 + 1
                    call cell_corners_grid(grid1, ic, jc, scx, scy)
                    call polygon_clip(scx, scy, tcx, tcy, ox, oy, no)
                    if (no >= 3) then
                        area = polygon_area(ox(1:no), oy(1:no))
                        if (area > 0.0_dp) then
                            nl = nl + 1
                            tsrc(nl) = cidx(cc)
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
        deallocate(emb, cidx, cd2, tsrc, tdst, tw)
        if (allocated(ox)) deallocate(ox, oy)
    end subroutine map_init_conservative

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
