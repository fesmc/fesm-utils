module root_finder

    use precision, only: wp

    implicit none

    private
    public :: solve_newton
    public :: solve_secant

contains

    subroutine solve_newton(x,n_iter,x_init,tol,n_max,f,fp,debug,converged)
        ! Estimate the zero of f(x) using Newton's method.
        ! Adapted from:
        ! https://faculty.washington.edu/rjl/classes/am583s2013/notes/fortran_newton.html

        implicit none

        real(wp), intent(OUT) :: x          ! Best guess of root
        integer,  intent(OUT) :: n_iter     ! Number of iterations to reach it
        real(wp), intent(IN)  :: x_init     ! Initial guess
        real(wp), intent(IN)  :: tol        ! Tolerance to convergence
        integer,  intent(IN)  :: n_max      ! Maximum iterations allowed
        real(wp), external    :: f          ! Function to find root of
        real(wp), external    :: fp         ! Derivative of Function
        logical,  intent(IN)  :: debug      ! Print iteration information?
        logical,  intent(OUT), optional :: converged  ! Did the method converge?

        ! Declare any local variables:
        real(wp) :: deltax, fx, fxprime
        integer  :: k
        logical  :: conv

        ! Save initial guess
        x    = x_init
        conv = .false.

        ! Newton iteration to find a zero of f(x)
        n_iter = 0

        do k = 1, n_max

            n_iter = n_iter + 1

            ! evaluate function and its derivative:
            fx      = f(x)
            fxprime = fp(x)

            if (debug) then
                write(*,*) n_iter, "x, f(x) = ", x, fx
            end if

            ! Converged if the residual is small enough
            if (abs(fx) < tol) then
                conv = .true.
                exit
            end if

            ! Guard against a zero derivative (flat point): the Newton step
            ! is undefined here, so stop rather than divide by zero.
            if (fxprime == 0.0_wp) then
                exit
            end if

            ! Compute Newton increment and update x:
            deltax = fx/fxprime
            x      = x - deltax

            ! Also converged if the step size is negligible. This catches
            ! roots where f is flat (e.g. multiple roots) and |f| stays
            ! above tol even though x has effectively stopped moving.
            if (abs(deltax) < tol) then
                conv = .true.
                exit
            end if

        end do

        if (present(converged)) converged = conv

        if (.not. conv) then
            write(*,*) "solve_newton:: Warning: no convergence."
        end if

        return

    end subroutine solve_newton

    subroutine solve_secant(x,n_iter,x_init,tol,n_max,f,debug,converged)
        ! Estimate the zero of f(x) using the Secant method.
        ! Adapted from:
        ! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/secant_f90.txt
        ! https://rosettacode.org/wiki/Roots_of_a_function#Fortran

        implicit none

        real(wp), intent(OUT) :: x          ! Best guess of root
        integer,  intent(OUT) :: n_iter     ! Number of iterations to reach it
        real(wp), intent(IN)  :: x_init     ! Initial guess
        real(wp), intent(IN)  :: tol        ! Tolerance to convergence
        integer,  intent(IN)  :: n_max      ! Maximum iterations allowed
        real(wp), external    :: f          ! Function to find root of
        logical,  intent(IN)  :: debug      ! Print iteration information?
        logical,  intent(OUT), optional :: converged  ! Did the method converge?

        ! Local variables
        integer  :: n
        real(wp) :: x1, x2
        real(wp) :: y1, y2
        real(wp) :: d
        logical  :: conv

        ! Set x to initial guess and a slightly different value
        x1 = x_init
        x2 = x_init*0.5_wp
        if (x_init == 0.0_wp) x2 = x_init + 0.1_wp

        x    = x2
        conv = .false.

        ! Evaluate f at the first point once; thereafter y1 is carried over
        ! from the previous iteration so f is evaluated only once per step.
        y1 = f(x1)

        n_iter = 0

        ! Start iterations
        do n = 1, n_max

            n_iter = n_iter + 1

            y2 = f(x2)

            if (debug) write(*,*) n_iter, x2, y2

            ! Converged if the residual is small enough
            if (abs(y2) < tol) then
                conv = .true.
                exit
            end if

            ! Guard against equal function values (flat secant): the update
            ! is undefined here, so stop rather than divide by zero.
            if (y2 == y1) then
                exit
            end if

            d = (x2 - x1) / (y2 - y1) * y2

            x1 = x2
            y1 = y2
            x2 = x2 - d

            ! Also converged if the step size is negligible.
            if (abs(d) < tol) then
                conv = .true.
                exit
            end if

        end do

        x = x2

        if (present(converged)) converged = conv

        if (.not. conv) then
            write(*,*) "solve_secant:: Warning: no convergence."
        end if

        return

    end subroutine solve_secant

end module root_finder
