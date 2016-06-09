module polar

use utility
use goal
use sensing

contains


subroutine getECPolar( i, N, p, cellCOM, ms, plrR, rSim, sigma, xcell, xCOM)
    implicit none
    integer,  intent(in) :: i, N
    real(b8), intent(in) :: ms, plrR
    real(b8), intent(inout), dimension(2) :: p
    real(b8), intent(in), dimension(:,:)  :: cellCOM
    integer,  intent(in), dimension(2)    :: rSim
    integer,  intent(in), dimension(:,:)  :: sigma, xcell
    real(b8), intent(in), dimension(2)    :: xCOM
    real(b8), dimension(2) :: q, qtmp
    integer,  dimension(N) :: nnL
    integer  :: j, P1, P2
    real(b8) :: msCOM, s

    q(:) = 0.0_b8

    ! get mean signal at cluster COM
    msCOM  = chemE( xCOM(1), xCOM(2))

    call getContactL( i, N, nnL, rSim, sigma, xcell)
    ! check if cell is enclosed
    P1 = int(perimCalc(rSim, sigma, xcell))
    P2 = sum(nnL)
    if( P1 /= P2 )then
        ! calculate repulsion vector
        do j = 1, N
            if( nnL(j) /= 0 )then
                qtmp = cellCOM(i,:) - cellCOM(j,:)
                qtmp = qtmp / sqrt( dot_product( qtmp, qtmp))
                q(:) = q(:) + real(nnL(j)) * qtmp
            endif
        enddo
        if( dot_product( q, q) /= 0.0 )then
            q = q / sqrt( dot_product( q, q))
        endif
        ! get signal value from distribution
        if( ms < 100.0 )then
            call poissonrand( ms, s)
        else
            s = normal(ms,sqrt(ms))
        endif
        if( s < 0.0 )then
            s = 0.0
        endif
        q = s * q / msCOM
    endif
    ! write(*,*) 'q =', q
    ! update polarization vector
    p = (1.0 - plrR) * p + (plrR*eps) * q

end subroutine getECPolar


! caclulate polarization vector using edge pixels of a cell
subroutine getMWPolar2( p, plrR, rSim, sigma, xcell)
    implicit none
    real(b8), intent(in) :: plrR
    integer,  intent(in), dimension(:,:)  :: sigma, xcell
    integer,  intent(in), dimension(2)    :: rSim
    real(b8), intent(inout), dimension(2) :: p
    real(b8), dimension(2) :: q, qtot
    integer,  dimension(2) :: nn
    integer :: edge, i, nl
    real(b8) :: ms, s, rx, ry

    qtot(:) = 0.0_b8
    nl = 1
    ! iterate over all pixels in a cell
    do while( xcell(nl,1) /= 0 )
        ! check whether pixel xcell(nl,:) is on an edge
        q(:) = 0.0_b8
        edge = 0
        do i = 1, 4
            call nnGet( i, nn, rSim, xcell(nl,:))
            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(xcell(nl,1),xcell(nl,2)) /= sigma(nn(1),nn(2)) )then
                ! ls is an edge
                q(:) = real(nn(:)) - real(xcell(nl,:)) + q(:)
                edge = edge + 1
            endif
        enddo
        ! normalize q vector and sample signal
        ms = 0.0
        s  = 0.0
        if( dot_product(q,q) /= 0.0 )then
            q = q / sqrt( dot_product(q,q) )
            ! write(*,*) 'pixel',nl,'q =',q
            rx = real(xcell(nl,1))
            ry = real(xcell(nl,2))
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
            ! write(*,*) 'ms =', ms,'s =',s,'q =', q
            q = s * q
            qtot = qtot + q
            ! write(*,*) 'qtot =',qtot
        endif

        nl = nl + 1
    enddo
    ! normalize qtot vector
    if( dot_product(qtot,qtot) /= 0.0 )then
        qtot = qtot / sqrt( dot_product(qtot,qtot) )
    endif
    ! write(*,*) 'qtot =', qtot
    ! update polarization vector
    p = (1.0 - plrR) * p + (plrR*eps) * qtot

end subroutine getMWPolar2


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
