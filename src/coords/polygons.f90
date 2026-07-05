module polygons
    ! Planar polygon geometry in double precision: area (shoelace),
    ! point-in-polygon (ray casting), and Sutherland-Hodgman clipping of a
    ! subject polygon against a convex clip polygon. These underpin the analytic
    ! conservative remapping (cell-overlap areas).

    use precision, only: dp

    implicit none
    private

    public :: polygon_area, point_in_polygon, polygon_clip
    public :: polygon_clip_sphere

contains

    subroutine polygon_clip_sphere(sv, cv, ov, no)
        ! Sutherland-Hodgman clipping on the unit sphere. Vertices are 3D unit
        ! vectors; polygon edges are great-circle arcs. The clip polygon cv must
        ! be convex and oriented counter-clockwise (interior on the +(a x b) side
        ! of each edge a->b). Output ov(3,:) has no vertices (0 if no overlap).
        real(dp), intent(in)  :: sv(:,:), cv(:,:)    ! (3, ns), (3, nc)
        real(dp), allocatable, intent(out) :: ov(:,:)
        integer,  intent(out) :: no

        real(dp), allocatable :: iv(:,:), wv(:,:)
        real(dp) :: a(3), b(3), p(3), q(3), x(3)
        integer  :: ns, nc, ni, nw, cap, e, k, kprev
        logical  :: pin, qin

        ns = size(sv,2); nc = size(cv,2)
        cap = ns + nc + 4
        allocate(iv(3,cap), wv(3,cap))
        ni = ns
        iv(:,1:ns) = sv

        do e = 1, nc
            a = cv(:,e); b = cv(:, mod(e,nc)+1)
            nw = 0
            if (ni == 0) exit
            do k = 1, ni
                kprev = k-1; if (kprev == 0) kprev = ni
                p = iv(:,kprev); q = iv(:,k)
                pin = (dot_product(q_cross(a,b), p) >= 0.0_dp)
                qin = (dot_product(q_cross(a,b), q) >= 0.0_dp)
                if (qin) then
                    if (.not. pin) then
                        call arc_intersect(p, q, a, b, x)
                        nw = nw+1; wv(:,nw) = x
                    end if
                    nw = nw+1; wv(:,nw) = q
                else if (pin) then
                    call arc_intersect(p, q, a, b, x)
                    nw = nw+1; wv(:,nw) = x
                end if
            end do
            ni = nw
            iv(:,1:nw) = wv(:,1:nw)
        end do

        no = ni
        allocate(ov(3, max(no,1)))
        if (no > 0) ov(:,1:no) = iv(:,1:no)
        deallocate(iv, wv)
    end subroutine polygon_clip_sphere

    pure function q_cross(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w(1) = u(2)*v(3) - u(3)*v(2)
        w(2) = u(3)*v(1) - u(1)*v(3)
        w(3) = u(1)*v(2) - u(2)*v(1)
    end function q_cross

    subroutine arc_intersect(s, e, a, b, x)
        ! Intersection of great-circle arc s->e with the great circle through a,b.
        real(dp), intent(in)  :: s(3), e(3), a(3), b(3)
        real(dp), intent(out) :: x(3)
        real(dp) :: nse(3), nab(3), dir(3), nrm
        nab = q_cross(a, b)
        nse = q_cross(s, e)
        dir = q_cross(nab, nse)
        nrm = sqrt(sum(dir**2))
        if (nrm < 1.0e-30_dp) then
            x = e; return
        end if
        dir = dir / nrm
        ! choose the intersection lying on the arc s->e (same hemisphere as its midpoint)
        if (dot_product(dir, s+e) < 0.0_dp) dir = -dir
        x = dir
    end subroutine arc_intersect

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
