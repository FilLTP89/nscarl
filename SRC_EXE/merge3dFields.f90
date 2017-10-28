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
    integer :: code, rk_, npr, comm; ! mpi communicator variables

    character(len=200) :: fld ! main working path
    integer :: nfm ! number of file to be read
    character(len=200), allocatable, dimension(:) :: fnm ! file names
    integer, dimension(:), allocatable :: dims,deltax

    integer :: i_,j_ ! counters

    ! attributes
    real(fpp), dimension(0:5) :: xLimBound
    integer(fpp), dimension(0:2) :: dimst
    real(fpp), dimension(0:2) :: xMinGlob, xMaxGlob

    ! partitioning
    real(fpp), dimension(0:2) :: xSplit
    
    comm = MPI_COMM_WORLD
    call init(comm)
    call init_hdf5()

    if(rk_ == 0) then
        write(*,*) 'Input file:',stdin
        read(*,'(I8)') nfm
        write(*,*) "nfiles = ", nfm
        allocate(fnm(0:nfm-1))
        allocate(dims(0:3*nfm-1))
        allocate(deltax(0:3*nfm-1))
        xLimBound(0:5) = 1e+20
        do i_=0,nfm-1
            read(*,'(A)') fnm(i_)
            fnm(i_) = "./"//trim(adjustL(fnm(i_)))  
            write(*,'(A)') 'FILE:',fnm(i_)
            call parse_mf_prop_nscarl(fnm(i_),xMinGlob,xMaxGlob,dimst)
            dims(3*i_:3*(i_+1)-1) = dimst
            write(*,*) 'DIMS:', dimst
            do j_=0,2
                deltax(3*i_+j_) = (xMaxGlob(j_)-xMinGlob(j_))/(dims(3*i_+j_)-1)
                xLimBound(0+j_) = min(xLimBound(0+j_),xMinGlob(j_))
                xLimBound(3+j_) = min(xLimBound(3+j_),xMaxGlob(j_))
            end do
        end do
        ! partitioning
        write(*,*) nint(npr/3.0)
        xSplit(0) = (xLimBound(3)-xLimBound(0))/max(nint(npr/3.0),1) 
        xSplit(1) = (xLimBound(4)-xLimBound(1))/max(nint(npr/3.0),1) 
        xSplit(2) = (xLimBound(5)-xLimBound(2))/(npr-2*max(nint(npr/3.0),1)) 
        write(*,*) 'Dims:',dims
        write(*,*) 'Min. Box Limits', xLimBound(0:2)
        write(*,*) 'Max. Box Limits', xLimBound(3:5)
        write(*,*) 'Delta X',deltax
        write(*,*) 'xSplit',xSplit
    end if
    
    !partitioning

    if(rk_ == 0) then
        deallocate(dims)
        deallocate(fnm)
    end if
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

