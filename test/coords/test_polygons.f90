program test_polygons
    ! Polygon area, point-in-polygon, and Sutherland-Hodgman clipping.

    use polygons
    use precision, only: dp

    implicit none

    real(dp) :: sx(4), sy(4), cx(4), cy(4)
    real(dp), allocatable :: ox(:), oy(:)
    integer :: no, fails
    real(dp) :: ar

    fails = 0

    ! Unit square area
    if (abs(polygon_area([0.0_dp,1.0_dp,1.0_dp,0.0_dp], &
                         [0.0_dp,0.0_dp,1.0_dp,1.0_dp]) - 1.0_dp) > 1.d-12) then
        write(*,*) "FAIL: unit square area"; fails=fails+1; end if

    ! Triangle area = 0.5
    if (abs(polygon_area([0.0_dp,1.0_dp,0.0_dp],[0.0_dp,0.0_dp,1.0_dp]) - 0.5_dp) > 1.d-12) then
        write(*,*) "FAIL: triangle area"; fails=fails+1; end if

    ! point-in-polygon
    if (.not. point_in_polygon(0.5_dp,0.5_dp, [0.0_dp,1.0_dp,1.0_dp,0.0_dp], &
                                              [0.0_dp,0.0_dp,1.0_dp,1.0_dp])) then
        write(*,*) "FAIL: point should be inside"; fails=fails+1; end if
    if (point_in_polygon(2.0_dp,2.0_dp, [0.0_dp,1.0_dp,1.0_dp,0.0_dp], &
                                        [0.0_dp,0.0_dp,1.0_dp,1.0_dp])) then
        write(*,*) "FAIL: point should be outside"; fails=fails+1; end if

    ! Overlap of [0,2]^2 and [1,3]^2 is [1,2]^2, area = 1
    sx = [0.0_dp,2.0_dp,2.0_dp,0.0_dp]; sy = [0.0_dp,0.0_dp,2.0_dp,2.0_dp]
    cx = [1.0_dp,3.0_dp,3.0_dp,1.0_dp]; cy = [1.0_dp,1.0_dp,3.0_dp,3.0_dp]
    call polygon_clip(sx, sy, cx, cy, ox, oy, no)
    ar = polygon_area(ox(1:no), oy(1:no))
    write(*,*) "overlap vertices =", no, " area =", ar, " (expect 1.0)"
    if (abs(ar - 1.0_dp) > 1.d-12) then; write(*,*) "FAIL: clip overlap area"; fails=fails+1; end if

    ! Subject fully inside clip -> area unchanged (0.25 square inside unit clip)
    call polygon_clip([0.25_dp,0.75_dp,0.75_dp,0.25_dp], [0.25_dp,0.25_dp,0.75_dp,0.75_dp], &
                      [0.0_dp,1.0_dp,1.0_dp,0.0_dp], [0.0_dp,0.0_dp,1.0_dp,1.0_dp], ox, oy, no)
    if (abs(polygon_area(ox(1:no),oy(1:no)) - 0.25_dp) > 1.d-12) then
        write(*,*) "FAIL: contained-polygon clip"; fails=fails+1; end if

    ! Disjoint -> no overlap
    call polygon_clip([0.0_dp,1.0_dp,1.0_dp,0.0_dp],[0.0_dp,0.0_dp,1.0_dp,1.0_dp], &
                      [5.0_dp,6.0_dp,6.0_dp,5.0_dp],[5.0_dp,5.0_dp,6.0_dp,6.0_dp], ox, oy, no)
    if (no /= 0) then; write(*,*) "FAIL: disjoint clip should be empty, got", no; fails=fails+1; end if

    if (fails > 0) then; write(*,*) "test_polygons FAILED with", fails; stop 1; end if
    write(*,*) "PASS: test_polygons"

end program test_polygons
