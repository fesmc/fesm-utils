program test_mapio
    ! Round-trip a weight_map through the SCRIP-superset NetCDF format.

    use coords

    implicit none

    type(weight_map_t) :: wm, wr
    integer :: fails

    fails = 0

    ! ---- MAP_DISTANCE: extensions must round-trip losslessly ----
    call weight_map_alloc(wm, MAP_DISTANCE, n_src=4, n_dst=2, n_links=3)
    wm%src = [1,2,3]; wm%dst = [1,1,2]
    wm%dist = [1.0_dp, 2.0_dp, 3.0_dp]
    wm%quadrant = [1,2,3]
    wm%xn = [0.1_dp, 0.2_dp, 0.3_dp]
    wm%yn = [0.4_dp, 0.5_dp, 0.6_dp]
    call weight_map_index(wm)
    call weight_map_write(wm, "test_wm_dist.nc")
    call weight_map_read(wr, "test_wm_dist.nc")

    if (wr%kind /= MAP_DISTANCE) then; write(*,*) "FAIL: kind not DISTANCE"; fails=fails+1; end if
    if (wr%n_src/=4 .or. wr%n_dst/=2 .or. wr%n_links/=3) then
        write(*,*) "FAIL: DISTANCE sizes"; fails=fails+1; end if
    if (any(wr%src/=wm%src) .or. any(wr%dst/=wm%dst)) then
        write(*,*) "FAIL: DISTANCE addresses"; fails=fails+1; end if
    if (maxval(abs(wr%dist-wm%dist))>1.d-12 .or. any(wr%quadrant/=wm%quadrant) .or. &
        maxval(abs(wr%xn-wm%xn))>1.d-12 .or. maxval(abs(wr%yn-wm%yn))>1.d-12) then
        write(*,*) "FAIL: DISTANCE payload"; fails=fails+1; end if
    call weight_map_free(wm); call weight_map_free(wr)

    ! ---- MAP_WEIGHT: standard SCRIP round-trip ----
    call weight_map_alloc(wm, MAP_WEIGHT, n_src=4, n_dst=2, n_links=3)
    wm%src = [1,2,3]; wm%dst = [1,1,2]; wm%w = [0.25_dp, 0.75_dp, 1.0_dp]
    call weight_map_index(wm)
    call weight_map_write(wm, "test_wm_wgt.nc")
    call weight_map_read(wr, "test_wm_wgt.nc")
    if (wr%kind /= MAP_WEIGHT) then; write(*,*) "FAIL: kind not WEIGHT"; fails=fails+1; end if
    if (any(wr%src/=wm%src) .or. any(wr%dst/=wm%dst) .or. maxval(abs(wr%w-wm%w))>1.d-12) then
        write(*,*) "FAIL: WEIGHT payload"; fails=fails+1; end if
    call weight_map_free(wm)

    ! ---- scrip_only on a distance map -> reload as frozen MAP_WEIGHT ----
    call weight_map_alloc(wm, MAP_DISTANCE, n_src=4, n_dst=2, n_links=3)
    wm%src=[1,2,3]; wm%dst=[1,1,2]; wm%dist=[1.0_dp,2.0_dp,3.0_dp]
    wm%quadrant=[1,2,3]; wm%xn=0.0_dp; wm%yn=0.0_dp
    call weight_map_index(wm)
    call weight_map_write(wm, "test_wm_so.nc", scrip_only=.true.)
    call weight_map_read(wr, "test_wm_so.nc")
    if (wr%kind /= MAP_WEIGHT) then
        write(*,*) "FAIL: scrip_only distance map should reload as WEIGHT"; fails=fails+1; end if
    ! dst1 baked shepard weights: (1/1, 1/4) normalized -> 0.8, 0.2
    if (abs(wr%w(1)-0.8_dp)>1.d-9 .or. abs(wr%w(2)-0.2_dp)>1.d-9) then
        write(*,*) "FAIL: baked shepard weights wrong", wr%w(1), wr%w(2); fails=fails+1; end if
    call weight_map_free(wm); call weight_map_free(wr)

    if (fails > 0) then
        write(*,*) "test_mapio FAILED with", fails, "errors"; stop 1
    end if
    write(*,*) "PASS: test_mapio"

end program test_mapio
