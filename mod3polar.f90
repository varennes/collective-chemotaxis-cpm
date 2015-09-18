module polar

use utility
use sensing

contains


! update the polarization vector
subroutine getPolar( p, plrR, comNew, comOld)
    implicit none
    real(b8), intent(in) :: plrR
    real(b8), intent(inout), dimension(2) :: p
    real(b8), intent(in),    dimension(2) :: comNew, comOld

    p = (1.0 - plrR) * p + COMnew - COMold

end subroutine getPolar


! update the polarization vector with sensing
subroutine getPolar2( p, plrR, eps, R0, Rk, comNew, comOld)
    implicit none
    real(b8), intent(in) :: plrR, eps, R0, Rk
    real(b8), intent(inout), dimension(2) :: p
    real(b8), intent(in),    dimension(2) :: comNew, comOld

    p = (1.0 - plrR) * p + eps * Rk/R0 * (COMnew - COMold)

end subroutine getPolar2


! update the polarization vector with sensing and antagonist mech, cooperative mech
subroutine getPolar3( p, plrR, q, R0, Rk, comNew, comOld)
    implicit none
    real(b8), intent(in) :: plrR, R0, Rk
    real(b8), intent(inout), dimension(2) :: p, q
    real(b8), intent(in),    dimension(2) :: comNew, comOld

    p = (1.0 - plrR) * p + (COMnew - COMold) + (eps*Rk/R0) * q

end subroutine getPolar3


! update the polarization vector with sensing and antagonist mech
subroutine getPolar4( p, plrR, q, R0, Rk)
    implicit none
    real(b8), intent(in) :: plrR, R0, Rk
    real(b8), intent(inout), dimension(2) :: p, q

    p = (1.0 - plrR) * p + (plrR*eps*Rk/R0) * q ! eps is now dimensionless

end subroutine getPolar4


! initialize the polarization vectors
subroutine itlPolar( N, plrP, p)
    implicit none
    integer, intent(in) :: N
    real(b8),    intent(in) :: plrP
    real(b8),    intent(out), dimension(:,:) :: p
    integer :: i
    real(b8) :: pi, r

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
    real(b8),    intent(out), dimension(:) :: com
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
real(b8) function getBias( aSig, bSig, plrP, p, x, xtmp)
    implicit none
    integer, intent(in) :: aSig, bSig
    real(b8),    intent(in) :: plrP
    real(b8),    intent(in), dimension(:,:)   :: p
    integer, intent(in), dimension(:,:,:) :: x, xtmp
    real(b8), dimension(2) :: comNew, comOld
    real(b8) :: sum1, sum2

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
                sum1 = sum1 / sqrt( dot_product( (comNew-comOld), (comNew-comOld)) )
            endif
        endif

        sum2 = 0.0
        if( bSig /= 0 )then
            call calcCellCOM(    x(bSig,:,:), comOld)
            call calcCellCOM( xtmp(bSig,:,:), comNew)

            sum2 = dot_product( (comNew-comOld), p(bSig,:))

            if( sum2 /= 0.0 )then
                sum2 = sum2 / sqrt( dot_product( p(bSig,:),p(bSig,:) ) )
                sum2 = sum2 / sqrt( dot_product( (comNew-comOld), (comNew-comOld)) )
            endif
        endif

        getBias = plrP * (sum1 + sum2)
    endif

end function getBias


! calculate bias due to polarization with sensing
real(b8) function getBias2( aSig, bSig, dXtMCS, p, x, xtmp)
    implicit none
    integer,  intent(in) :: aSig, bSig
    real(b8), intent(in), dimension(:)     :: dXtMCS
    real(b8), intent(in), dimension(:,:)   :: p
    integer,  intent(in), dimension(:,:,:) :: x, xtmp
    real(b8), dimension(2) :: comNew, comOld
    real(b8) :: sum1, sum2

    getBias2 = 0.0

    sum1 = 0.0
    if( aSig /= 0 .AND. dXtMCS(aSig) > 1e-10 )then
        call calcCellCOM(    x(aSig,:,:), comOld)
        call calcCellCOM( xtmp(aSig,:,:), comNew)

        sum1 = dot_product( (comNew-comOld), p(aSig,:))

        if( sum1 /= 0.0 )then
            sum1 = sum1 / dXtMCS(aSig)
            sum1 = sum1 / sqrt( dot_product( (comNew-comOld), (comNew-comOld)) )
        endif
    endif

    sum2 = 0.0
    if( bSig /= 0 .AND. dXtMCS(bSig) > 1e-10 )then
        call calcCellCOM(    x(bSig,:,:), comOld)
        call calcCellCOM( xtmp(bSig,:,:), comNew)

        sum2 = dot_product( (comNew-comOld), p(bSig,:))

        if( sum2 /= 0.0 )then
            sum2 = sum2 / dXtMCS(bSig)
            sum2 = sum2 / sqrt( dot_product( (comNew-comOld), (comNew-comOld)) )
        endif
    endif

    getBias2 = sum1 + sum2

end function getBias2


end module
