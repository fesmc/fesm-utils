program test

    use precision
    
    !use mapping_scrip
    use gaussian_quadrature

    type(gq2D_class) :: gq2D
    type(gq3D_class) :: gq3D

    real(wp) :: var_qp(4)
    real(wp) :: var(2,2)
    real(wp) :: dx
    real(wp) :: dy 
    integer  :: i, j 
    integer  :: im1, ip1, jm1, jp1

    logical, parameter :: test_2D = .TRUE.
    logical, parameter :: test_3D = .FALSE.
    
    if (test_2D) then
        call gq2D_init(gq2D,verbose=.TRUE.)

        var(1,1) = 2.0
        var(1,2) = 2.0
        var(2,1) = 4.0
        var(2,2) = 4.0

        dx = 1.0
        dy = 1.0 
        
        i = 2
        j = 2
        im1 = i-1
        ip1 = i
        jm1 = j-1
        jp1 = j

        call gq2D_to_nodes(gq2D, var_qp, var, dx, dy, "ab",i,j,im1,ip1,jm1,jp1)

        write(*,*) "var_qp: ", var_qp
    
    end if

    write(*,*)
    write(*,*) "----"
    write(*,*)

    if (test_3D) then
        call gq3D_init(gq3D,verbose=.TRUE.)
    end if

end program test
