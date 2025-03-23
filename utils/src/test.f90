program test

    !use mapping_scrip
    use gaussian_quadrature

    type(gq2D_class) :: gq2D
    type(gq3D_class) :: gq3D

    real(wp) :: var_qp(4)
    real(wp) :: var(2,2)
    real(wp) :: dx
    real(wp) :: dy 

if (.FALSE.) then
    call gq2D_init(gq2D)

    write(*,*)
    write(*,*) "----"
    write(*,*)

    call gq3D_init(gq3D)
end if

    ! ==== CISM style ====

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
    
end program test
