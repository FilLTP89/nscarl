module h5field_parse

    implicit none

contains
    subroutine read_dset_3d_real(parent, name, data, ibase)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        real(fpp), dimension(:,:,:), allocatable, intent(out) :: data
        integer, intent(in), optional :: ibase
        !
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(3) :: dims, maxdims
        integer :: hdferr, i0, i1,i2
        i0 = 1
        i1 = 1
        i2 = 1
        if (present(ibase)) then
            i0 = ibase
            i1 = ibase
            i2 = ibase
        end if

        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_3d_real : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_3d_real : h5dgetspace KO"
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr .ne. 3) stop "read_dset_3d_real : h5sgetdim KO "
        allocate(data(i0:dims(1)+i0-1, i1:dims(2)+i1-1, i2:dims(3)+i2-1))
        call h5dread_f(dset_id, H5T_REAL, data, dims, hdferr)
        if (hdferr .ne. 0) then
            write(*,*) "read_dset_3d_real : h5dread KO", dims
            stop 1
        endif
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_3d_real : h5dclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_3d_real : h5sclose KO"
    end subroutine read_dset_3d_real

    subroutine read_subset_3d_real(parent, name, imin, imax, data)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        integer, intent(in), dimension(0:2) :: imin, imax
        real(fpp), dimension(:,:,:), allocatable, intent(out) :: data
        !
        integer(HID_T) :: dset_id, space_id, memspace_id
        integer(HSIZE_T), dimension(3) :: start, count
        integer :: hdferr

        start = imin
        count = imax-imin+[1,1,1]
        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5dgetspace KO"
        call H5Sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, start, count, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5sgetdim KO "
        allocate(data(imin(0):imax(0), imin(1):imax(1), imin(2):imax(2)))
        call H5Screate_simple_f(3, count, memspace_id, hdferr)
        call h5dread_f(dset_id, H5T_REAL, data, count, hdferr, memspace_id, space_id)
        if (hdferr .ne. 0) then
            write(*,*) "read_subset_3d_real : h5dread KO", imin, imax
            stop 1
        end if
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5dclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5sclose KO"
        call h5sclose_f(memspace_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5sclose KO"
    end subroutine read_subset_3d_real

    subroutine read_attr_real_vec(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        real(fpp), dimension(:), intent(out) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims, max_dims

        dims(1) = size(data)
        call h5aopen_f(file_id, attr_name, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        call h5sget_simple_extent_dims_f(space_id, dims, max_dims, hdferr) !Get the size of the attribute
        call h5aread_f(attr_id, H5T_REAL, data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine read_attr_real_vec
    subroutine read_attr_real(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        real(fpp), intent(out) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        dims(1) = 1
        call h5aopen_f(dset, attr, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        ! TODO: check h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        call h5aread_f(attr_id, H5T_REAL, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine read_attr_real

    subroutine write_attr_real(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        real(fpp), intent(in) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims

        dims(1) = 1
        !write(*,*) "save_attr_real: ", attr, value
        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
        call h5acreate_f(dset, attr, H5T_NATIVE_DOUBLE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_REAL, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_attr_real

    subroutine write_dataset_d3(parent, name, arr)
        use HDF5
        implicit none
        integer(HID_T), intent(in) :: parent
        character(len=*), intent(in) :: name
        real(fpp), dimension(:,:,:), intent(in) :: arr
        !
        integer(HSIZE_T), dimension(3) ::  dims
        integer(HID_T) :: dset_id
        integer :: hdferr
        dims(1) = size(arr,1)
        dims(2) = size(arr,2)
        dims(3) = size(arr,3)
        call create_dset_3d_i8(parent, name, H5T_IEEE_F64LE, dims, dset_id)
        call h5dwrite_f(dset_id, H5T_REAL, arr, dims, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_d3 : h5dwrite KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_d3 : h5dclose KO"
    end subroutine write_dataset_d3

    subroutine create_dset_3d_i8(parent, name, dtype, dims, dset_id)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent, dtype
        integer(HID_T), INTENT(OUT) :: dset_id
        integer(HSIZE_T), INTENT(INOUT), dimension(3) :: dims
        integer(HSIZE_T), dimension(3) :: chunk, maxdims
        integer(HID_T) :: space_id, prop_id
        integer :: hdferr

        maxdims = dims
        chunk(1) = min(dims(1), 256_HSIZE_T)
        chunk(2) = min(dims(2), 1024_HSIZE_T)
        if (dims(3)==H5S_UNLIMITED_F) then
            chunk(3) = 64
            dims(3) = 0
        else
            chunk(3) = max(1_HSIZE_T, min(dims(3), int(256*1024/(chunk(1)*chunk(2)),HSIZE_T)))
        endif
        call h5screate_simple_f(3, dims, space_id, hdferr, maxdims)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5screate KO"
        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pcreate KO"
        if ((dims(1)*dims(2)*dims(3)).gt.128 .or. dims(3)==H5S_UNLIMITED_F) then
            call h5pset_chunk_f(prop_id, 3, chunk, hdferr)
            if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_chunk KO"
            if (dtype/=H5T_IEEE_F32LE .and. dtype/=H5T_IEEE_F64LE) then
                ! Les donnees en float donnent un taux de compression bas pour un
                ! cout de calcul eleve
                call h5pset_deflate_f(prop_id, 5, hdferr)
                if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_deflate KO"
                call h5pset_shuffle_f(prop_id, hdferr)
                if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_shuffle KO"
            endif
        end if
        call h5dcreate_f(parent, name, dtype, space_id, dset_id, hdferr, prop_id)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5dcreate KO"
        call h5pclose_f(prop_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5sclose KO"
    end subroutine create_dset_3d_i8
end module h5field_parse
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=4 et tw=80 smartindent :
