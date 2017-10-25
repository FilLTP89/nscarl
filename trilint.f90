module trilint

    implicit none

contains

    subroutine interpolate_elem_field(Tdomain, specel, mat, pf, field)
        type(domain), intent(inout) :: Tdomain
        type(element), intent(inout) :: specel
        type(subdomain), intent(in) :: mat
        type(PropertyField), intent(in) :: pf
        !
        real(fpp), dimension(0:mat%ngll-1,0:mat%ngll-1,0:mat%ngll-1) :: field
        real(fpp), dimension(0:2) :: xx ! node coord
        real(fpp), dimension(0:2) :: aa   ! node interpolation coeffs
        real(fpp) :: pos, val
        integer :: i,j,k     ! gll index
        integer :: idef      ! global coord index
        integer :: n
        integer, dimension(0:2) :: ii  ! index of x,y,z inside material domain
        logical :: pml

        pml = is_pml(mat)
        do k = 0,mat%ngll-1
            do j = 0,mat%ngll-1
                do i = 0,mat%ngll-1

                    idef = specel%Iglobnum(i,j,k)
                    do n=0,2
                        xx(n) = Tdomain%GlobCoord(n,idef)
                        ! Traitement PML : on echantillonne au bord de la PML uniquement
                        if (pml) then
                            if (mat%pml_width(n)>0) then
                                if (xx(n)>mat%pml_pos(n)) xx(n)=mat%pml_pos(n)
                            end if
                            if (mat%pml_width(n)<0) then
                                if (xx(n)<mat%pml_pos(n)) xx(n)=mat%pml_pos(n)
                            end if
                        end if
                        pos = (xx(n)-pf%MinBound(n))/pf%step(n)
                        ii(n) = floor(pos)
                        if (ii(n)  < pf%imin(n)) ii(n) = pf%imin(n)
                        if (ii(n) >= pf%imax(n)) ii(n) = pf%imax(n)-1
                        aa(n) = pos-ii(n)
                    end do
                    ! trilinear interpolation
                    val =       (1.-aa(0))*(1.-aa(1))*(1.-aa(2))*pf%var(ii(0)  ,ii(1)  ,ii(2)  )
                    val = val + (   aa(0))*(1.-aa(1))*(1.-aa(2))*pf%var(ii(0)+1,ii(1)  ,ii(2)  )
                    val = val + (1.-aa(0))*(   aa(1))*(1.-aa(2))*pf%var(ii(0)  ,ii(1)+1,ii(2)  )
                    val = val + (   aa(0))*(   aa(1))*(1.-aa(2))*pf%var(ii(0)+1,ii(1)+1,ii(2)  )
                    val = val + (1.-aa(0))*(1.-aa(1))*(   aa(2))*pf%var(ii(0)  ,ii(1)  ,ii(2)+1)
                    val = val + (   aa(0))*(1.-aa(1))*(   aa(2))*pf%var(ii(0)+1,ii(1)  ,ii(2)+1)
                    val = val + (1.-aa(0))*(   aa(1))*(   aa(2))*pf%var(ii(0)  ,ii(1)+1,ii(2)+1)
                    val = val + (   aa(0))*(   aa(1))*(   aa(2))*pf%var(ii(0)+1,ii(1)+1,ii(2)+1)
                    field(i,j,k) = val
                end do
            end do
        end do
    end subroutine interpolate_elem_field

end module trilint
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
