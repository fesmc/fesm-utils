program test_roots
    ! Exercise the root_finder module (solve_newton, solve_secant) against
    ! several functions with known roots, plus degenerate cases that must
    ! fail gracefully (zero derivative / flat secant) rather than crash.

    use precision, only: wp
    use root_finder

    implicit none

    integer,  parameter :: n_max = 100
    real(wp), parameter :: tol   = real(1e-6,wp)   ! solver tolerance on |f|/|dx|
    real(wp), parameter :: xtol  = real(1e-3,wp)   ! accepted error on the root

    integer :: nfail

    nfail = 0

    write(*,*) "==================================================="
    write(*,*) " root_finder test suite"
    write(*,*) "==================================================="
    write(*,"(a)") ""
    write(*,"(a4,2x,a8,2x,a18,2x,a12,2x,a12,2x,a5,2x,a4)") &
        "meth","func","x (computed)","x (true)","err","iter","ok"

    ! --- Simple roots: both methods should find them accurately ---

    ! x^2 - 2 = 0  ->  x = sqrt(2)
    call check_newton("x^2-2 ", f_sq,  fp_sq,  1.5_wp, sqrt(2.0_wp))
    call check_secant("x^2-2 ", f_sq,          2.0_wp, sqrt(2.0_wp))

    ! cos(x) - x = 0  ->  Dottie number ~ 0.7390851
    call check_newton("cos-x ", f_cx,  fp_cx,  0.5_wp, 0.7390851_wp)
    call check_secant("cos-x ", f_cx,          0.5_wp, 0.7390851_wp)

    ! exp(x) - 2 = 0  ->  x = ln(2)
    call check_newton("exp-2 ", f_ex,  fp_ex,  1.0_wp, log(2.0_wp))
    call check_secant("exp-2 ", f_ex,          1.0_wp, log(2.0_wp))

    ! (x+1)^5 = 0  ->  multiple root at x = -1 (flat near root). Convergence
    ! is only linear here and (x+1)^5 underflows near the root in single
    ! precision, so both methods stop well before machine accuracy. Assert a
    ! deliberately loose tolerance. The secant start avoids x_init*0.5 = -1
    ! landing exactly on the root, so the iteration is genuinely exercised.
    call check_newton("(x+1)^5", f_p5, fp_p5, -2.0_wp, -1.0_wp, 0.1_wp)
    call check_secant("(x+1)^5", f_p5,        -1.6_wp, -1.0_wp, 0.1_wp)

    ! --- Degenerate cases: must report non-convergence, not crash ---

    ! f(x) = 1 has no root; derivative is 0 everywhere.
    call check_nonconv_newton("no-root", f_one, fp_zero, 0.0_wp)
    call check_nonconv_secant ("no-root", f_one,          0.0_wp)

    write(*,"(a)") ""
    write(*,"(a)") "==================================================="
    if (nfail == 0) then
        write(*,"(a)") " All root_finder tests PASSED."
    else
        write(*,"(a,i0,a)") " root_finder tests FAILED: ", nfail, " case(s)."
    end if
    write(*,"(a)") "==================================================="

    if (nfail > 0) stop 1

contains

    subroutine check_newton(name,f,fp,x_init,x_true,xtol_in)
        character(len=*), intent(IN) :: name
        real(wp), external :: f, fp
        real(wp), intent(IN) :: x_init, x_true
        real(wp), intent(IN), optional :: xtol_in
        real(wp) :: x, err, xt
        integer  :: n_iter
        logical  :: converged, ok

        xt = xtol
        if (present(xtol_in)) xt = xtol_in

        call solve_newton(x,n_iter,x_init,tol,n_max,f,fp,.false.,converged)
        err = abs(x - x_true)
        ok  = converged .and. err < xt
        call report("newt",name,x,x_true,err,n_iter,ok)
        if (.not. ok) nfail = nfail + 1
    end subroutine check_newton

    subroutine check_secant(name,f,x_init,x_true,xtol_in)
        character(len=*), intent(IN) :: name
        real(wp), external :: f
        real(wp), intent(IN) :: x_init, x_true
        real(wp), intent(IN), optional :: xtol_in
        real(wp) :: x, err, xt
        integer  :: n_iter
        logical  :: converged, ok

        xt = xtol
        if (present(xtol_in)) xt = xtol_in

        call solve_secant(x,n_iter,x_init,tol,n_max,f,.false.,converged)
        err = abs(x - x_true)
        ok  = converged .and. err < xt
        call report("sec ",name,x,x_true,err,n_iter,ok)
        if (.not. ok) nfail = nfail + 1
    end subroutine check_secant

    subroutine check_nonconv_newton(name,f,fp,x_init)
        ! Expect converged == .false. and a finite (non-NaN) result.
        character(len=*), intent(IN) :: name
        real(wp), external :: f, fp
        real(wp), intent(IN) :: x_init
        real(wp) :: x
        integer  :: n_iter
        logical  :: converged, ok

        call solve_newton(x,n_iter,x_init,tol,n_max,f,fp,.false.,converged)
        ok = (.not. converged) .and. (x == x)   ! x==x is false for NaN
        call report_flag("newt",name,"expect no-conv",ok)
        if (.not. ok) nfail = nfail + 1
    end subroutine check_nonconv_newton

    subroutine check_nonconv_secant(name,f,x_init)
        character(len=*), intent(IN) :: name
        real(wp), external :: f
        real(wp), intent(IN) :: x_init
        real(wp) :: x
        integer  :: n_iter
        logical  :: converged, ok

        call solve_secant(x,n_iter,x_init,tol,n_max,f,.false.,converged)
        ok = (.not. converged) .and. (x == x)
        call report_flag("sec ",name,"expect no-conv",ok)
        if (.not. ok) nfail = nfail + 1
    end subroutine check_nonconv_secant

    subroutine report(meth,name,x,x_true,err,n_iter,ok)
        character(len=*), intent(IN) :: meth, name
        real(wp), intent(IN) :: x, x_true, err
        integer,  intent(IN) :: n_iter
        logical,  intent(IN) :: ok
        character(len=4) :: tag
        tag = "PASS"
        if (.not. ok) tag = "FAIL"
        write(*,"(a4,2x,a8,2x,f18.10,2x,f12.7,2x,es12.3,2x,i5,2x,a4)") &
            meth, name, x, x_true, err, n_iter, tag
    end subroutine report

    subroutine report_flag(meth,name,msg,ok)
        character(len=*), intent(IN) :: meth, name, msg
        logical,          intent(IN) :: ok
        character(len=4) :: tag
        tag = "PASS"
        if (.not. ok) tag = "FAIL"
        write(*,"(a4,2x,a8,2x,a,2x,a4)") meth, name, msg, tag
    end subroutine report_flag

    ! --- Test functions and their derivatives ---

    real(wp) function f_sq(x)
        real(wp) :: x
        f_sq = x*x - 2.0_wp
    end function f_sq
    real(wp) function fp_sq(x)
        real(wp) :: x
        fp_sq = 2.0_wp*x
    end function fp_sq

    real(wp) function f_cx(x)
        real(wp) :: x
        f_cx = cos(x) - x
    end function f_cx
    real(wp) function fp_cx(x)
        real(wp) :: x
        fp_cx = -sin(x) - 1.0_wp
    end function fp_cx

    real(wp) function f_ex(x)
        real(wp) :: x
        f_ex = exp(x) - 2.0_wp
    end function f_ex
    real(wp) function fp_ex(x)
        real(wp) :: x
        fp_ex = exp(x)
    end function fp_ex

    real(wp) function f_p5(x)
        real(wp) :: x
        f_p5 = (x + 1.0_wp)**5
    end function f_p5
    real(wp) function fp_p5(x)
        real(wp) :: x
        fp_p5 = 5.0_wp*(x + 1.0_wp)**4
    end function fp_p5

    real(wp) function f_one(x)
        real(wp) :: x
        f_one = 1.0_wp
    end function f_one
    real(wp) function fp_zero(x)
        real(wp) :: x
        fp_zero = 0.0_wp
    end function fp_zero

end program test_roots
