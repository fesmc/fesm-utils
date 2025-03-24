program test

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

if (.TRUE.) then
    call gq2D_init(gq2D)

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
    
    write(*,*)
    write(*,*) "----"
    write(*,*)

    !call gq3D_init(gq3D)
end if

    ! ==== CISM style ====
if (.FALSE.) then
    write(*,*)
    write(*,*) "----------------------------------------"
    write(*,*)
    write(*,*) "CISM style..."
    write(*,*)

    call gaussian_quadrature_init()

    var(1,1) = 2.0
    var(1,2) = 2.0
    var(2,1) = 4.0
    var(2,2) = 4.0

    dx = 1.0
    dy = 1.0 

    call gaussian_quadrature_2D_to_nodes(var_qp,var,dx,dy,i=2,j=2,grid_type="ab")

    write(*,*) "var_qp: ", var_qp
end if

end program test
