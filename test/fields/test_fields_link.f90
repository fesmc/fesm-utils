program test_fields_link
    ! fesm_fields_check_link, accepted cases: equivalent units spelled
    ! differently, different horizontal grids, and the four allowed rows of
    ! the time-semantics table (INSTANT->INSTANT, MEAN->INSTANT, MEAN->MEAN,
    ! ACCUM->ACCUM). Each rejected case (units, rank, nz, cons and the five
    ! forbidden time-semantics pairs) is an abort program in abort/.

    use fesm_fields

    implicit none

    type(field_state_t) :: src, dst
    integer :: hs(4), hd(4), hsmb, hdsmb

    call fields_init(src, "test:link_src")
    call fields_init(dst, "test:link_dst")

    ! Compatible pair: same quantity in two spellings, on two grids.
    call fields_add(src, "smb", "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, hsmb,  grid="gs", dims=[4,3])
    call fields_add(dst, "smb", "mm w.e./yr",  FIELD_INSTANT, FIELD_EXTENSIVE, hdsmb, grid="gd", dims=[8,6])

    ! The allowed rows of the time-semantics table.
    call fields_add(src, "a", "K", FIELD_INSTANT, FIELD_INTENSIVE, hs(1), grid="gs", dims=[4,3,12], zdim="month")
    call fields_add(dst, "a", "K", FIELD_INSTANT, FIELD_INTENSIVE, hd(1), grid="gd", dims=[8,6,12], zdim="month")
    call fields_add(src, "b", "K", FIELD_MEAN,    FIELD_INTENSIVE, hs(2))
    call fields_add(dst, "b", "K", FIELD_INSTANT, FIELD_INTENSIVE, hd(2))
    call fields_add(src, "c", "K", FIELD_MEAN,    FIELD_INTENSIVE, hs(3))
    call fields_add(dst, "c", "K", FIELD_MEAN,    FIELD_INTENSIVE, hd(3))
    call fields_add(src, "d", "kg m-2", FIELD_ACCUM, FIELD_EXTENSIVE, hs(4), grid="gs", dims=[4,3])
    call fields_add(dst, "d", "mm we",  FIELD_ACCUM, FIELD_EXTENSIVE, hd(4), grid="gd", dims=[8,6])

    ! Links are checked on the declaration; locking is not required.
    call fields_check_link(src, hsmb, dst, hdsmb)
    write(*,"(a)") "  ok    MEAN kg m-2 yr-1 (4x3) -> INSTANT mm w.e./yr (8x6)"
    call fields_check_link(src, hs(1), dst, hd(1))
    write(*,"(a)") "  ok    INSTANT -> INSTANT, rank 3"
    call fields_check_link(src, hs(2), dst, hd(2))
    write(*,"(a)") "  ok    MEAN -> INSTANT"
    call fields_check_link(src, hs(3), dst, hd(3))
    write(*,"(a)") "  ok    MEAN -> MEAN"
    call fields_check_link(src, hs(4), dst, hd(4))
    write(*,"(a)") "  ok    ACCUM -> ACCUM"

    call fields_lock(src)
    call fields_lock(dst)
    call fields_check_link(src, hsmb, dst, hdsmb)
    write(*,"(a)") "  ok    same link after lock"

    write(*,*) "PASS: test_fields_link"

end program test_fields_link
