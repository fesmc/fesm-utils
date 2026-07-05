program test_proj
    ! Projection check — define a polar-stereographic point set (Bamber corners)
    ! and inverse-project x/y to lon/lat via the coords/oblimap machinery.

    use coords

    implicit none

    type(points_class) :: pts
    integer :: i

    call points_init(pts, name="Bamber-corners", mtype="polar_stereographic", units="kilometers", &
                     lambda=-39.d0, phi=90.d0, alpha=19.d0, lon180=.TRUE., &
                     x=[-800.d0,-800.d0,700.d0,700.d0], y=[-3400.d0,-600.d0,-600.d0,-3400.d0])

    write(*,"(a)") "        x         y    =>        lat       lon"
    do i = 1, size(pts%x,1)
        write(*,"(2f10.3,a,2f10.3)") pts%x(i), pts%y(i), " => ", pts%lat(i), pts%lon(i)
    end do

    ! Sanity: all corners should land in the northern hemisphere over Greenland.
    if (any(pts%lat < 55.d0) .or. any(pts%lat > 85.d0)) then
        write(*,*) "FAIL: projected latitudes outside expected Greenland range"
        stop 1
    end if
    write(*,*) "PASS: test_proj"

end program test_proj
