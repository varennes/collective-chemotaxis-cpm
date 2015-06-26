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
        p(i,1) = plrP * cos( 2.0 * pi * r )
        p(i,2) = plrP * sin( 2.0 * pi * r )
    enddo

end subroutine itlPolar


! calculate center of mass of a single cell
subroutine calcCellCOM( xcell, com)
    ! x = array of all the cell lattice sites
    ! com = center of mass of a single cell
    implicit none
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


! calculate bias due to polarization
real function getBias( aSig, bSig, plrP, p, x, xtmp)
    implicit none
    integer, intent(in) :: aSig, bSig
    real,    intent(in) :: plrP
    real,    intent(in), dimension(:,:)   :: p
    integer, intent(in), dimension(:,:,:) :: x, xtmp
    real, dimension(2) :: comNew, comOld
    real :: sum1, sum2

    ! check if plrP = 0.0
    if( plrP == 0.0 )then
        getBias = 0.0
    else
        sum1 = 0.0
        if( aSig /= 0 )then
            call calcCellCOM(    x(aSig,:,:), comOld)
            call calcCellCOM( xtmp(aSig,:,:), comNew)

            sum1 = dot_product( (comNew-comOld), p(aSig,:))

            if( sum1 /= 0.0 )then
                sum1 = sum1 / sqrt( dot_product( p(aSig,:),p(aSig,:) ) )
            endif
        endif

        sum2 = 0.0
        if( bSig /= 0 )then
            call calcCellCOM(    x(bSig,:,:), comOld)
            call calcCellCOM( xtmp(bSig,:,:), comNew)

            sum2 = dot_product( (comNew-comOld), p(bSig,:))

            if( sum2 /= 0.0 )then
                sum2 = sum2 / sqrt( dot_product( p(bSig,:),p(bSig,:) ) )
            endif
        endif

        getBias = plrP * (sum1 + sum2)
    endif

end function getBias


end module
