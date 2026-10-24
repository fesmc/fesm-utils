program nml_keys
    ! List every key of the namelist files given on the command line, one
    ! line `file group key` per key, groups and keys in file order
    ! (nml_group_names, nml_group_keys). The key list read by
    ! test_nml_compat (test/nml/nml_compat.sh).

    use nml, only : nml_group_names, nml_group_keys, nml_key_len, nml_group_len

    implicit none

    character(len=2048) :: fname
    character(len=nml_group_len), allocatable :: groups(:)
    character(len=nml_key_len),   allocatable :: keys(:)
    integer :: i, ig, ik

    do i = 1, command_argument_count()
        call get_command_argument(i, fname)
        groups = nml_group_names(trim(fname))
        do ig = 1, size(groups)
            keys = nml_group_keys(trim(fname), trim(groups(ig)))
            do ik = 1, size(keys)
                write(*,"(a)") trim(fname)//" "//trim(groups(ig))//" "//trim(keys(ik))
            end do
        end do
    end do

end program nml_keys
