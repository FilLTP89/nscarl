program merge3dFields

    use mpi
    use hdf5
    use constants
    use sem_hdf5
    use matStructure

#ifdef f2003
    use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
        stdout=>output_unit, &
        stderr=>error_unit
#else
#define stdin 5
#define stdout 6
#define stderr 0
#endif

    implicit none

    !INPUTS
    integer :: code, rk_, npr, comm ! mpi communicator variables

    character(len=200) :: fld ! main working path
    integer :: nfm ! number of file to be read
    character(len=200), allocatable, dimension(:) :: fnm,prop ! file names
    integer, dimension(:), allocatable :: dims

    ! counters 
    integer :: i_,j_,k_
    integer :: nf,np
    ! attributes
    real(fpp), dimension(0:5) :: xLimBound
    integer(fpp), dimension(0:2) :: dimst, npt
    real(fpp), dimension(0:2) :: xMinGlob, xMaxGlob

    ! partitioning
    real(fpp), dimension(0:2) :: dxSplit
    real(fpp), dimension(:), allocatable :: xSplit, deltax  
    ! data samples
    integer(fpp), dimension(:,:), allocatable :: idxg
    real(fpp), allocatable, dimension(:,:,:) :: datasamples

    comm = MPI_COMM_WORLD
    call init(comm)
    call init_hdf5()

    if(rk_ == 0) then
        write(*,*) 'Input file:',stdin
        read(*,'(I8)') nfm
        write(*,*) "nfiles = ", nfm
    end if    
    call MPI_BCAST(nfm, 1, MPI_INTEGER, 0, comm, code) 
    allocate(fnm(0:nfm-1))
    allocate(prop(0:nfm-1))
    if(rk_ == 0) then
        nf=0
        np=0
        do i_=0,nfm-1
            read(*,'(A)') fnm(i_)
            fnm(i_) = "./"//trim(adjustL(fnm(i_)))  
            nf=nf+len(fnm(i_))
            write(*,*) 'FILE:',fnm(i_)
            read(*,*) prop(i_) 
            prop(i_) = trim(adjustL(prop(i_)))  
            np=np+len(prop(i_))
            write(*,*) 'PROP:',prop(i_)
        end do
        close(stdin,status='keep')
    end if    
    call MPI_BCAST(nf, 1, MPI_INTEGER, 0, comm, code) 
    call MPI_BCAST(np, 1, MPI_INTEGER, 0, comm, code) 
    call MPI_BCAST(fnm, nf, MPI_CHARACTER, 0, comm, code) 
    call MPI_BCAST(prop, np, MPI_CHARACTER, 0, comm, code) 
    allocate(dims(0:3*nfm-1))
    allocate(deltax(0:3*nfm-1))
    xLimBound(0:5) = 1e+20
    do i_=0,nfm-1
        call parse_mf_prop_nscarl(fnm(i_),xMinGlob,xMaxGlob,dimst)
        dims(3*i_:3*(i_+1)-1) = dimst
        if (rk_.eq.0) write(*,*) 'DIMS:', dimst
        do j_=0,2
            deltax(3*i_+j_) = (xMaxGlob(j_)-xMinGlob(j_))/(dims(3*i_+j_)-1)
            xLimBound(0+j_) = min(xLimBound(0+j_),xMinGlob(j_))
            xLimBound(3+j_) = min(xLimBound(3+j_),xMaxGlob(j_))
        end do
    end do
    ! partitioning
    npt(0) = ceiling(npr**(1.0/3.0))
    npt(1) = npt(0) 
    npt(2) = npr-product(npt(0:1))
    if (rk_==0) write(*,*) 'npt',npt
    allocate(xSplit(0:npr+2))
    k_=0
    do i_=0,2
        dxSplit(i_) = (xLimBound(3+i_)-xLimBound(i_))/npt(i_)
        do j_=0,npt(i_)
            xSplit(k_+j_) = xLimBound(i_)+dxSplit(i_)*j_
        end do
        k_=k_+npt(i_)+1
    end do
    allocate(idxg(0:2,0:npr-1))
    if (rk_.eq.0)  call create_global_index(npr,npt,3,idxg)
    if (rk_.eq.0) write(*,*) idxg
    call MPI_BARRIER(comm,code)
    deallocate(dims)
    deallocate(fnm)
    call finalize()

    contains

        subroutine init(comm_local)
            implicit none
            !INPUT
            integer, intent(in) :: comm_local

            call MPI_INIT(code)
            call MPI_COMM_RANK(comm_local, rk_, code)
            call MPI_COMM_SIZE(comm_local, npr, code)
        end subroutine init

        subroutine finalize()
            implicit none
            call MPI_FINALIZE(code)
        end subroutine finalize

end program merge3dFields
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
