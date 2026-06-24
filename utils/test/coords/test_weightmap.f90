program test_weightmap
    ! Unit-test the weight_map store + apply kernel on hand-built maps.

    use weight_map
    use precision, only: dp

    implicit none

    type(weight_map_t) :: wm
    real(dp) :: var1(3), var2w(2), var2d(1)
    logical  :: m2(2)
    integer  :: fails
    real(dp), parameter :: miss = -9999.0_dp

    fails = 0
    var1 = [10.0_dp, 20.0_dp, 30.0_dp]

    ! ---- MAP_WEIGHT: dst1 <- src1(w=1),src2(w=3) ; dst2 <- src3(w=1) ----
    call weight_map_alloc(wm, MAP_WEIGHT, n_src=3, n_dst=2, n_links=3)
    wm%dst = [1,1,2]; wm%src = [1,2,3]; wm%w = [1.0_dp,3.0_dp,1.0_dp]
    call weight_map_index(wm)

    call weight_map_apply(wm, var1, var2w, "mean", missing_value=miss, mask2=m2)
    call check("WEIGHT mean dst1", var2w(1), 17.5_dp)      ! (10+60)/4
    call check("WEIGHT mean dst2", var2w(2), 30.0_dp)
    if (.not. (m2(1) .and. m2(2))) then; write(*,*) "FAIL: mask2"; fails=fails+1; end if

    call weight_map_apply(wm, var1, var2w, "count", missing_value=miss)
    call check("WEIGHT count dst1 (mode=20)", var2w(1), 20.0_dp)
    call weight_map_free(wm)

    ! ---- MAP_DISTANCE: dst1 <- src1(d=1,q=1), src2(d=2,q=2) ----
    call weight_map_alloc(wm, MAP_DISTANCE, n_src=3, n_dst=1, n_links=2)
    wm%dst = [1,1]; wm%src = [1,2]; wm%dist = [1.0_dp,2.0_dp]
    wm%quadrant = [1,2]; wm%xn = 0.0_dp; wm%yn = 0.0_dp
    call weight_map_index(wm)

    call weight_map_apply(wm, var1, var2d, "nn", missing_value=miss)
    call check("DISTANCE nn (nearest=10)", var2d(1), 10.0_dp)

    call weight_map_apply(wm, var1, var2d, "shepard", missing_value=miss)
    call check("DISTANCE shepard", var2d(1), 12.0_dp)      ! (10*1 + 20*0.25)/1.25

    call weight_map_apply(wm, var1, var2d, "quadrant", missing_value=miss)
    call check("DISTANCE quadrant", var2d(1), 12.0_dp)

    ! missing source -> nn falls through to next valid neighbor
    var1(1) = miss
    call weight_map_apply(wm, var1, var2d, "nn", missing_value=miss)
    call check("DISTANCE nn with masked nearest (->20)", var2d(1), 20.0_dp)
    call weight_map_free(wm)

    if (fails > 0) then
        write(*,*) "test_weightmap FAILED with", fails, "errors"; stop 1
    end if
    write(*,*) "PASS: test_weightmap"

contains

    subroutine check(label, got, expect)
        character(len=*), intent(in) :: label
        real(dp),         intent(in) :: got, expect
        if (abs(got - expect) > 1.0e-10_dp) then
            write(*,*) "FAIL: ", label, " got", got, " expected", expect
            fails = fails + 1
        else
            write(*,*) "ok:   ", label, " =", got
        end if
    end subroutine check

end program test_weightmap
