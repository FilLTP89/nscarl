module matStructure
    use sem_hdf5
    use constants
    implicit none

    type LMC_properties

        ! variables d'écrouissage kinematic et isotrope de Lamaitre et Chaboche
        real(fpp) :: sigma_yld   ! first yielding limit
        real(fpp) :: C_kin       ! variable for kinematic hardening
        real(fpp) :: kapa_kin    ! variable for kinematic hardening
        real(fpp) :: b_iso       ! variable for isotropic hardening
        real(fpp) :: Rinf_iso    ! variable for isotropic hardening

    end type LMC_properties
    !
    type HYP_properties

        ! hypoelastique model
        real(fpp) :: a0,a1,a2 ! loading 
        real(fpp) :: b0,b1,b2 ! unloading

    end type HYP_properties
    !
    type nl_properties
        type(LMC_properties) :: LMC_prop
        type(HYP_properties) :: HYP_prop
    end type nl_properties
    !

    type PropertyField
        character(len=1024) :: propFilePath
        character(len=100) :: propName ! name of property and HDF5 group
        real(fpp), dimension(0:2) :: MinBound, MaxBound, step
        ! XXX: need to handle fields of different sizes...
        integer, dimension(0:2) :: NN ! dimension of the grid for this property
        integer, dimension(0:2) :: imin, imax
        real(fpp), dimension(:,:,:), allocatable :: var
        type(nl_properties) :: nl_prop
    end type PropertyField

    type Subdomain
        integer          :: dom ! The computation domain SOLID/FLUID/SPML/FPML
        integer          :: material_definition
        integer          :: deftype
        logical          :: present ! true if an element with this mat exists on this cpu
        !! Numerotation gll
        integer :: NGLL

        !! Definition materiau solide, isotrope
        real(fpp) :: Pspeed, Sspeed, Ddensity
        real(fpp) :: DLambda, DMu
        real(fpp) :: DE, DNu
        real(fpp) :: DKappa
        real(kind=8) :: dt
        !! Definition materiau solide anisotrope
        ! TODO

        !! NONLINEAR LEMAITRE-CHABOCHE
        ! real(fpp) :: DSyld,DCkin,DKkin
        real(fpp) :: DRinf,DBiso
        real(fpp) :: DNlkp
        !! HYPOELASTIC
        real(fpp) :: da0,da1,da2,db0,db1,db2 
        
        !! ATTENUATION
        real(fpp) :: Qmu, Qpression
        !! PML
        real(fpp), dimension(0:2) :: pml_pos, pml_width
        integer :: npow
        real(fpp) :: Apow

        !! Boundaries for material initialisation from file
        real(fpp), dimension(0:2) :: MinBound_Loc, MaxBound_Loc

        ! three variables on a 3D grid used for intializing from fields in files
        ! depending on material_definition we can have
        ! Vp(v1) Vs(v2) Rho(v3)
        ! Lambda(v1) Mu(v2) Rho(v3) ...
        type(PropertyField), dimension(3) :: pf

    end type Subdomain

    private :: gcoord, gindex

contains
    ! etant donne une grille de n points regulierement espaces entre xmin et xmax
    ! revoie grille(i)
    function gcoord(i, n, xmin, xmax) result(x)
        integer, intent(in) :: i,n
        real(fpp), intent(in) :: xmin, xmax
        real(fpp) :: x

        x = xmin + i*(xmax-xmin)/(n-1)
    end function gcoord

    ! soit une grille de n points reguliers entre xmin et xmax, renvoie
    ! i tel que x(i)<=x<x(i+1)   x(0) = xmin x(n-1)=xmax
    function gindex(x, n, xmin, xmax) result(i)
        integer, intent(in) :: n
        real(fpp), intent(in) :: x, xmin, xmax
        integer :: i

        i = floor( ((n-1)*(x-xmin))/(xmax-xmin) )
    end function gindex

    function prop_check_var(pid, pname, gid)
        use hdf5
        integer(HID_T), intent(in) :: pid
        integer(HID_T), intent(out) :: gid
        character(len=100) :: pname
        !
        logical :: ok, prop_check_var
        integer :: hdferr
        !
        prop_check_var = .false.
        call H5Lexists_f(pid, trim(adjustl(pname)), ok, hdferr)
        if (ok) then
            call H5Gopen_f(pid, trim(adjustl(pname)), gid, hdferr)
            if (hdferr == 0) prop_check_var = .true.
        endif
    end function prop_check_var

    subroutine parse_mf_prop_nscarl(fnm,xMinGlob,xMaxGlob,dims)

        character(len=200), intent(in) :: fnm ! file names
        real(fpp), dimension(0:2), intent(inout) :: xMinGlob, xMaxGlob
        integer(HSIZE_T), dimension(0:2), intent(inout) :: dims
        integer(HSIZE_T), dimension(:), allocatable :: dimst
        integer ::hdferr ! hdf5 error variable
        integer(HID_T) :: fid ! hdf5 file identifier

        call h5fopen_f(fnm, H5F_ACC_RDONLY_F, fid, hdferr)
        call read_attr_real_vec(fid, 'xMinGlob',xMinGlob(0:2))
        call read_attr_real_vec(fid, 'xMaxGlob',xMaxGlob(0:2))
        call read_dims(fid, "samples", dimst)
        dims=dimst 
        call h5fclose_f(fid, hdferr)
    end subroutine

    subroutine init_prop_file_field(mat, pf)
        use hdf5
        use sem_hdf5
        type(subdomain), intent(inout) :: mat
        type(PropertyField), intent(inout) :: pf
        !
        integer :: k, hdferr
        integer(HID_T) :: file_id, grp_id
        integer(HSIZE_T), dimension(:), allocatable :: dims
        logical :: subgrp

        if (.not. mat%present) return
        call init_hdf5()
        call h5fopen_f(trim(pf%propFilePath), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
        if(hdferr /= 0) then
            write(*,*) "Could not open file:", trim(pf%propFilePath)
            stop 1
        end if
        subgrp = prop_check_var(file_id, pf%propName, grp_id)
        if (.not. subgrp) then
            grp_id = file_id
        end if
        call read_attr_real_vec(grp_id, "xMinGlob", pf%MinBound)
        call read_attr_real_vec(grp_id, "xMaxGlob", pf%MaxBound)
        call read_dims(grp_id, "samples", dims)
        pf%NN = dims
        ! On va calculer les indices i0,i1 j0,j1,k0,k1 tels que
        ! i0 plus grand entier tel que x(i0)<MinBound_loc(0), 
        ! i1 plus petit entier tel que x(i1)>MaxBound_loc(0), etc... 

        do k = 0,2
            pf%imin(k) = gindex(mat%MinBound_Loc(k), pf%NN(k), pf%MinBound(k), pf%MaxBound(k))
            pf%imax(k) = gindex(mat%MaxBound_Loc(k), pf%NN(k), pf%MinBound(k), pf%MaxBound(k))+1
            pf%step(k) = (pf%MaxBound(k)-pf%MinBound(k))/(pf%NN(k)-1)
            if (pf%imin(k)<0) pf%imin(k) = 0
            if (pf%imax(k)<pf%imin(k)) pf%imax(k) = pf%imin(k)
            if (pf%imax(k)>=pf%NN(k)) pf%imax(k) = pf%NN(k)-1
            if (pf%imin(k)>pf%imax(k)) pf%imin(k) = pf%imax(k)
        end do
        
        call read_subset_3d_real(grp_id, "samples", pf%imin, pf%imax, pf%var)
        if (subgrp) call H5Gclose_f(grp_id, hdferr)
        call H5Fclose_f(file_id, hdferr)
    end subroutine init_prop_file_field
end module matStructure
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
