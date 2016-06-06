module polar

use utility
use sensing

contains


    ! update polarization vector for Many Wrongs (MW) mechanism
    subroutine getMWPolar( p, plrR, cellCOM, xcell)
        implicit none
        real(b8), intent(in) :: plrR
        integer,  intent(in),  dimension(:,:) :: xcell
        real(b8), intent(in),  dimension(:)   :: cellCOM
        real(b8), intent(inout), dimension(2) :: p
        real(b8), dimension(2) :: q, qpixel
        integer :: i, j, nl
        real(b8) :: ms, s, rx, ry, qmag

        q(:) = 0.0
        qmag = 1.0
        call occupyCount( nl, xcell(:,:) )
        ! iterate over all pixels within one cell
        do i = 1, nl
            ms = 0.0
            s  = 0.0
            qpixel(:) = 0.0

            rx = real(xcell(i,1))
            ry = real(xcell(i,2))
            ! get mean signal at that point
            ms  = chemE( rx, ry)
            ! get signal value from distribution
            if( ms < 100.0 )then
                call poissonrand( ms, s)
            else
                s = normal(ms,sqrt(ms))
            endif
            if( s < 0.0 )then
                s = 0.0
            endif

            ! calculate pixel vector
            qpixel(1) = rx - cellCOM(1)
            qpixel(2) = ry - cellCOM(2)

            q = q + s * qpixel
            ! write(*,*) ' s = ',s,' q =',q
            ! write(*,*) ' qpx = ',qpixel
        enddo
        ! make q a unit vector
        qmag = sqrt( dot_product( q, q) )
        ! write(*,*) ' qmag =', qmag
        q = q / qmag
        ! write(*,*) ' q =', q, sqrt( dot_product( q, q) )

        p = (1.0 - plrR) * p + (plrR*eps) * q

    end subroutine getMWPolar


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
    ! write(*,*) 'getBias2 function', aSig, bsig

    sum1 = 0.0
    if( aSig /= 0 .AND. dXtMCS(aSig) > 1e-10 )then
        call calcCellCOM(    x(aSig,:,:), comOld)
        call calcCellCOM( xtmp(aSig,:,:), comNew)

        sum1 = dot_product( (comNew-comOld), p(aSig,:))
        ! write(*,*) 'p =', p(aSig,:)

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
        write(*,*) 'p =', p(aSig,:)

        if( sum2 /= 0.0 )then
            sum2 = sum2 / dXtMCS(bSig)
            sum2 = sum2 / sqrt( dot_product( (comNew-comOld), (comNew-comOld)) )
        endif
    endif

    getBias2 = sum1 + sum2

end function getBias2


end module
