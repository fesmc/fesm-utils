module polygons
    ! Planar polygon geometry in double precision: area (shoelace),
    ! point-in-polygon (ray casting), and Sutherland-Hodgman clipping of a
    ! subject polygon against a convex clip polygon. These underpin the analytic
    ! conservative remapping (cell-overlap areas).

    use precision, only: dp

    implicit none
    private

    public :: polygon_area, point_in_polygon, polygon_clip

contains

    pure function polygon_area(x, y) result(area)
        ! Unsigned area of the polygon with vertices (x,y) via the shoelace formula.
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: area
        integer  :: n, i, j
        n = size(x)
        area = 0.0_dp
        if (n < 3) return
        j = n
        do i = 1, n
            area = area + (x(j) + x(i)) * (y(j) - y(i))
            j = i
        end do
        area = 0.5_dp * abs(area)
    end function polygon_area

    pure function signed_area(x, y) result(area)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: area
        integer  :: n, i, j
        n = size(x)
        area = 0.0_dp
        j = n
        do i = 1, n
            area = area + (x(j)*y(i) - x(i)*y(j))
            j = i
        end do
        area = 0.5_dp * area
    end function signed_area

    pure function point_in_polygon(px, py, x, y) result(inside)
        ! Ray-casting test: is (px,py) inside polygon (x,y)?
        real(dp), intent(in) :: px, py, x(:), y(:)
        logical :: inside
        integer :: n, i, j
        n = size(x)
        inside = .false.
        j = n
        do i = 1, n
            if (((y(i) > py) .neqv. (y(j) > py))) then
                if (px < (x(j)-x(i)) * (py-y(i)) / (y(j)-y(i)) + x(i)) then
                    inside = .not. inside
                end if
            end if
            j = i
        end do
    end function point_in_polygon

    subroutine polygon_clip(sx, sy, cx, cy, ox, oy, no)
        ! Sutherland-Hodgman: clip subject polygon (sx,sy) against the CONVEX
        ! clip polygon (cx,cy). Output polygon (ox,oy) has no vertices (0 if no
        ! overlap). ox,oy must be allocatable; they are sized here.
        real(dp), intent(in)  :: sx(:), sy(:), cx(:), cy(:)
        real(dp), allocatable, intent(out) :: ox(:), oy(:)
        integer,  intent(out) :: no

        real(dp), allocatable :: ix(:), iy(:)   ! input list
        real(dp), allocatable :: wx(:), wy(:)   ! working (output) list
        real(dp) :: ax, ay, bx, by              ! current clip edge a->b
        real(dp) :: ccx(size(cx)), ccy(size(cy))
        integer  :: nc, ns, ni, nw, cap, e, k, kprev
        real(dp) :: px, py, qx, qy
        logical  :: pin, qin

        ns = size(sx); nc = size(cx)

        ! Ensure clip polygon is counter-clockwise (interior on the left)
        if (signed_area(cx, cy) < 0.0_dp) then
            do k = 1, nc
                ccx(k) = cx(nc-k+1); ccy(k) = cy(nc-k+1)
            end do
        else
            ccx = cx; ccy = cy
        end if

        cap = ns + nc + 4
        allocate(ix(cap), iy(cap), wx(cap), wy(cap))

        ! start with the subject polygon as the input list
        ni = ns
        ix(1:ns) = sx; iy(1:ns) = sy

        do e = 1, nc
            ax = ccx(e); ay = ccy(e)
            bx = ccx(mod(e, nc) + 1); by = ccy(mod(e, nc) + 1)
            nw = 0
            if (ni == 0) exit
            do k = 1, ni
                kprev = k - 1; if (kprev == 0) kprev = ni
                px = ix(kprev); py = iy(kprev)
                qx = ix(k);     qy = iy(k)
                pin = inside_edge(px, py, ax, ay, bx, by)
                qin = inside_edge(qx, qy, ax, ay, bx, by)
                if (qin) then
                    if (.not. pin) then
                        nw = nw + 1
                        call seg_line_intersect(px, py, qx, qy, ax, ay, bx, by, wx(nw), wy(nw))
                    end if
                    nw = nw + 1
                    wx(nw) = qx; wy(nw) = qy
                else
                    if (pin) then
                        nw = nw + 1
                        call seg_line_intersect(px, py, qx, qy, ax, ay, bx, by, wx(nw), wy(nw))
                    end if
                end if
            end do
            ni = nw
            ix(1:nw) = wx(1:nw); iy(1:nw) = wy(1:nw)
        end do

        no = ni
        allocate(ox(max(no,1)), oy(max(no,1)))
        if (no > 0) then
            ox(1:no) = ix(1:no); oy(1:no) = iy(1:no)
        end if
        deallocate(ix, iy, wx, wy)
    end subroutine polygon_clip

    pure logical function inside_edge(px, py, ax, ay, bx, by) result(inside)
        ! Is point (px,py) on the interior (left) side of directed edge a->b?
        real(dp), intent(in) :: px, py, ax, ay, bx, by
        inside = ((bx-ax)*(py-ay) - (by-ay)*(px-ax)) >= 0.0_dp
    end function inside_edge

    subroutine seg_line_intersect(px, py, qx, qy, ax, ay, bx, by, ix, iy)
        ! Intersection of segment p->q with the infinite line through a,b.
        real(dp), intent(in)  :: px, py, qx, qy, ax, ay, bx, by
        real(dp), intent(out) :: ix, iy
        real(dp) :: r1, r2, denom, t
        r1 = (bx-ax)*(py-ay) - (by-ay)*(px-ax)
        r2 = (bx-ax)*(qy-ay) - (by-ay)*(qx-ax)
        denom = r1 - r2
        if (abs(denom) < 1.0e-30_dp) then
            ix = qx; iy = qy
        else
            t = r1 / denom
            ix = px + t*(qx-px)
            iy = py + t*(qy-py)
        end if
    end subroutine seg_line_intersect

end module polygons
