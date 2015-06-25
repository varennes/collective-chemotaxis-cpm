module polar

use utility

contains


! update the polarization vector
subroutine getPolar( p, plrR, comNew, comOld)
    implicit none
    real, intent(in) :: plrR
    real, intent(inout), dimension(2) :: p
    real, intent(in),    dimension(2) :: comNew, comOld

    p = (1.0 - plrR) * p + COMnew - COMold

end subroutine getPolar


! initialize the polarization vectors
subroutine itlPolar( N, plrP, p)
    implicit none
    integer, intent(in) :: N
    real,    intent(in) :: plrP
    real,    intent(out), dimension(:,:) :: p
    integer :: i
    real :: pi, r

    pi = 3.14159265359
    p  = 0.0

    do i = 1, N
        call random_number(r)
        p(1) = plrP * cos( 2.0 * pi * r )
        p(1) = plrP * sin( 2.0 * pi * r )
    enddo

end subroutine itlPolar


! calculate center of mass of a single cell
subroutine calcCellCOM( xcell, com)
    ! x = array of all the cell lattice sites
    ! com = center of mass of a single cell
    implicit none
    integer, intent(in) :: N
    integer, intent(in),  dimension(:,:) :: xcell
    real,    intent(out), dimension(:) :: com
    integer :: i, j, nl

    com = 0.0

    call occupyCount( nl, xcell(:,:))

    do i = 1, nl
        do j = 1, 2
            com(j) = com(j) + real(xcell(i,j))
        enddo
    enddo

    com = com / real(nl)

end subroutine calcCellCOM


end module
