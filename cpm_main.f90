program cpmcollab

use utility
use goal
use polar
use sensing
use simpleconnect
use wrtout

! allocate variables
implicit none
integer :: i, j, k, ne, ne0, nf, nl, check, aSig, bSig
integer :: tcount, tELEM, tMCS, tmax
integer :: nRun, runTotal

integer :: A0, N, x1, x2, P1, P2, clusterP, clusterA
integer, dimension(2) :: r0, rCell, rSim

integer, allocatable :: sigma(:,:), sigmaTmp(:,:)
integer, allocatable :: x(:,:,:), xTmp(:,:,:)
integer, allocatable :: edge(:,:), filled(:,:), node(:), nnL(:,:)
integer, dimension(2) :: a, b, nn

real(b8), allocatable :: dXtMCS(:), p(:,:), q(:,:), cellCOM(:,:), cellCOMold(:,:), ptmp(:,:), deltaCOM(:)
real(b8) :: plrP, plrR, P0

real(b8), allocatable :: firstpass(:),  MSD(:), MSDrun(:), xCOM(:,:)
real(b8) :: d, df, dreset, prob, r, uNew, uOld, w
real(b8) :: neMean, neMeanRun
real(b8) :: t0, tf, treset

real(b8) :: rx1, rx2, speciesR0
real(b8), allocatable  :: signal(:), meanSignal(:), speciesR(:), speciesX(:), speciesY(:), meanY(:)
real(b8), allocatable  :: etaY(:), gNN(:,:), M(:,:)
real(b8), dimension(2) :: qtmp

call cpu_time(t0)

! initialize parameters, variables
open(unit=11,file='input.txt',status='old',action='read')
read(11,*) r0(1), r0(2)
read(11,*) rCell(1), rCell(2)
read(11,*) x1, x2   ! rSim(1) = x1+x2, simulation size (length)
read(11,*) rSim(2)  ! rSim(2) simulation size (width)
read(11,*) tmax     ! max number of time steps simulation will run
read(11,*) runTotal ! total number of runs that will be simulated
read(11,*) df       ! threshold distance for recording first-passage time

rSim(1) = x1 + x2
N       = rCell(1) * rCell(2) ! total number of cells
A0      = r0(1) * r0(2)       ! relaxed cell size
P0      = 3.6*sqrt( real(A0)) ! relaxed cell perimeter
! speciesR0 = g * sqrt( real(N) * real(A0)**3.0 )
dreset = 20.0

! initialize parameters for polarization
open(unit=12,file='polarInput.txt',status='old',action='read')
read(12,*) plrP ! initalized polarization vector magnitude
read(12,*) plrR ! polarization vector decay rate

write(*,*) '   N =',N
write(*,*) '  A0 =',A0, ' P0 =', P0
write(*,*) 'rSim =',rSim
write(*,*) '  df =',df
write(*,*) 'plrP =',plrP,' plrR =',plrR
write(*,*)

call init_random_seed()

allocate(    sigma( rSim(1) + 2, rSim(2) + 2) )
allocate( sigmaTmp( rSim(1) + 2, rSim(2) + 2) )
allocate(     x( N, 4*A0, 2) )
allocate(  xTmp( N, 4*A0, 2) )
allocate(  xCOM( tmax, 2) )

allocate( edge( (rSim(1)+2)*(rSim(2)+2), 2) )
allocate( filled( 4*A0, 2) )
allocate( node(2) )

allocate(    MSD( tmax) )
allocate( MSDrun( tmax) )
allocate( firstpass(runTotal) )

allocate( p(N,2) )
allocate( q(N,2) )

allocate( ptmp(N,2) )

allocate( cellCOM(N,2) )
allocate( cellCOMold(N,2) )
allocate( dXtMCS(N) )
allocate( deltaCOM(N*(N-1)/2) )
allocate( signal(N) )
allocate( meanSignal(N) )
allocate( speciesR(N) )
allocate( speciesX(N) )
allocate( meanY(N) )
allocate( speciesY(N) )
allocate(  etaY(N) )
allocate( gNN(N,N) )
allocate(   M(N,N) )
allocate( nnL(N,N) )

MSDrun    = 0.0
neMeanrun = 0.0
firstpass = 0.0

do i = 2, rSim(1)+1
    rx1 = real(i)
    rx2 = real(rSim(2))
    write(122,*) chemE( rx1, rx2)
enddo

do nRun = 1, runTotal

write(*,*) '  nRun =',nRun

tcount = 1
tELEM  = 1 ! elementary time step
tMCS   = 1 ! MCS step, 1 MCS step = L^2 ELEM steps
treset = 0.0

edge   = 0
node   = 0
filled = 0

x    = 0
xTmp = 0
dXtMCS   = 0.0
sigma    = 0
sigmaTmp = 0

d    = 0.0
xCOM = 0.0
MSD  = 0.0

a = 1
b = 1

! initialize sigma
rSim(1) = x1
call itlSigma( r0, rCell, rSim, sigma)
! call itlSigmaRandom( N, r0, rCell, rSim, sigma)
! call itlSigmaSuper( N, r0, rCell, rSim, sigma)

! initialize edge
call itlEdge( edge, ne, rSim, sigma)
ne0 = ne
neMean = real(ne)

! initialize x array
rSim(1) = x1
call makeX( N, rSim, sigma, x)

! initialize xCOM and cellCOM
! initialize polarization
call calcXCOM( N, x, xCOM(1,:))
do i = 1, N
    call calcCellCOM( x(i,:,:),  cellCOM(i,:))
    ! call getMWPolar( p(i,:), plrR, cellCOM(i,:), x(i,:,:))
    call getMWPolar2( p, plrR, rSim, sigma, x(i,:,:))
enddo
cellCOMold = cellCOM

write(*,*) ' beta =',beta
write(*,*) ' eps  =',eps
write(*,*)

! calculate initial energy
rSim(1) = x1 + x2
uOld = goalEval1( A0, P0, N, rSim, sigma, x)
uNew = 0.0
! write(*,*) ' uNew=', uNew, 'uOld=', uOld

! write outputs
! call wrtSigma( rSim, sigma, tELEM)
! call wrtEdge( edge, rSim, sigma, tELEM)
! call wrtEdgeArray( edge, tELEM)
! call wrtU( 0.0, uOld, 0.0, 0.0, tELEM)
! call wrtXR( N, x, speciesR, tELEM)
call wrtPolar( N, p, tELEM)! call wrtX( N, x, tELEM)
! call wrtX( N, x, tELEM)
! do i = 1, N
!     write(155,*) cellCOM(i,:), tELEM - 1
! enddo


do while( tMCS < tmax )
    tELEM = tELEM + 1

    call countEdge( edge, ne)

    call picka( a, edge, ne)

    call pickb( a, b, rSim)

    if( b(1) == 1 .OR. b(2) == 1 )then
        ! don't attempt to copy sigma value
    elseif( b(1) == (rSim(1) + 2) .OR. b(2) == (rSim(2) + 2) )then
        ! don't attempt to copy sigma value

    elseif( sigma(b(1),b(2)) /= sigma(a(1),a(2)) )then
        ! make and update sigmaTmp
        sigmaTmp = sigma
        sigmaTmp(b(1),b(2)) = sigma(a(1),a(2))

        call makeX( N, rSim, sigmaTmp, xTmp)

        ! check if the new configuration is simply connected
        check = 0
        do i = 1, N
            filled = 0
            node = xTmp(i,1,:)

            call floodFill( node, filled, rSim, xTmp(i,:,:))

            call occupyCount( nl, xTmp(i,:,:))
            call occupyCount( nf, filled)

            if( nl == nf .AND. nl /= 0 )then
                check = check + 1
            endif
        enddo

        if( check == N )then
            ! the new configuration is simply connected
            ! calculate probability of execution
            aSig = sigma(a(1),a(2))
            bSig = sigma(b(1),b(2))

            w = getBias2( aSig, bSig, dXtMCS, p, x, xtmp)

            ! write(*,*) '  w =',w,'a =',aSig,'b =',bSig,'dx =',dXtMCS(aSig),dXtMCS(bSig)

            uNew = goalEval1( A0, P0, N, rSim, sigmaTmp, xTmp)
            prob = probEval( uNew, uOld, w)

            ! write(*,*)
            ! write(*,*) ' uNew=', uNew, 'uOld=', uOld, ' w=', w

            call random_number(r)

            if( r < prob )then
                ! execute
                sigma = sigmaTmp
                x     = xTmp
                uOld  = uNew

                ! check whether a and b nn's are/still an edge
                do i = 1, 4
                    call nnGet( i, nn, rSim, a)

                    if( nn(1) /= 0 )then
                        call updateEdge( nn, edge, ne, rSim, sigma)
                    endif

                    call nnGet( i, nn, rSim, b)

                    if( nn(1) /= 0 )then
                        call updateEdge( nn, edge, ne, rSim, sigma)
                    endif
                enddo

            endif ! execution  if statement
        endif ! simple connect if statement

    endif

    ! check if simulation has reached a MCS time step
    if( mod( tELEM - 1, ne0) == 0 )then
        tMCS   = tMCS   + 1
        tcount = tcount + 1

        neMean = neMean + real(ne)

        k = 0
        deltaCOM = 0.0
        do i = 1, N
            call calcCellCOM( x(i,:,:),  cellCOM(i,:))
            ! calculate intercell distances
            do j = i+1, N
                k = k + 1
                deltaCOM(k) = sqrt( dot_product( cellCOM(i,:)-cellCOM(j,:), cellCOM(i,:)-cellCOM(j,:) ) )
            enddo
            dXtMCS(i) = sqrt( dot_product( cellCOM(i,:)-cellCOMold(i,:), cellCOM(i,:)-cellCOMold(i,:) ) )
        enddo

        ! update polarization vector
        do i = 1, N
            call calcCellCOM( x(i,:,:),  cellCOM(i,:))
            ! write(*,*) 'cellCOM =', cellCOM(i,:)
            ! call getMWPolar( p(i,:), plrR, cellCOM(i,:), x(i,:,:))
            call getMWPolar2( p, plrR, rSim, sigma, x(i,:,:))
        enddo

        cellCOMold = cellCOM

        ! calculate COM of the whole group
        call calcXCOM( N, x, xCOM(tMCS,:))
        ! calculate MSD
        MSD(tMCS) = calcMSD( xCOM(tMCS,:), xCOM(1,:))

        ! write outputs
        if( mod( tMCS-1, 10) == 0)then

            ! output cluster area and perimeter
            ! if( treset /= 0 )then
            !     clusterA = 0
            !     clusterP = 0
            !     P1 = 0
            !     P2 = 0
            !     do i = 1, N
            !         P1 = int(perimCalc(rSim, sigma, x(i,:,:))) ! whole cell perimeter
            !         P2 = sum(nnL(i,:)) ! total cell-cell contact
            !         if( P2 /= 0 )then
            !             call occupyCount( nl, x(i,:,:) )
            !             clusterP = P1 - P2 + clusterP
            !             clusterA = nl + clusterA
            !         endif
            !     enddo
            !     write(190,*) clusterA, clusterP, tMCS-1
            ! endif

            ! write(161,'(I7)', advance='no') tMCS-1 ! write species X
            ! do i = 1, N
            !     write(161,'(F9.2)', advance='no') speciesX(i)
            ! enddo
            ! write(161,*) ''
            ! write(162,'(I7)', advance='no') tMCS-1 ! write species Y
            ! do i = 1, N
            !     write(162,'(F9.2)', advance='no') speciesY(i)
            ! enddo
            ! write(162,*) ''

            ! call wrtSigma( rSim, sigma, tMCS)
            ! write(150,*) xCOM(tMCS,:), tMCS
            ! call wrtXR( N, x, speciesR, tMCS)
            call wrtPolar( N, p, tMCS)
            call wrtX( N, x, tMCS)
            do i = 1, N
                write(155,*) cellCOM(i,:), tMCS - 1
            enddo
        endif

        ! calculate d
        d = calcD( xCOM(tMCS,1), xCOM(1,1))
        if( treset == 0.0 )then
            if( d >= dreset .AND. d < (dreset + 2.0) )then
                treset = tMCS - 1
            endif
        endif
        if( d >= df )then
            firstpass(nRun) = tMCS - 1 - treset
            tMCS = tmax
        endif

    endif

enddo ! end while loop

MSDrun = MSDrun + MSD

neMean    = neMean / real(tcount)

write(108,*) firstpass(nRun) * real(ne0)/neMean

neMeanRun = neMeanRun + neMean

write(*,*) '  FPT: treset =', treset

enddo ! end run loop

MSDrun = MSDrun / real(runTotal)
! do i = 1, tmax
!     write(107,*) MSDrun(i), i-1
! enddo

neMeanRun = neMeanRun / real(runTotal)

do i = 1, runTotal
    write(109,*) firstpass(i) * real(ne0)/neMeanRun
enddo


close(11)

call cpu_time(tf)
write(*,*)
write(*,*) 'Run time =',tf-t0
write(*,*) 'Initial edge area =',ne0
write(*,*) 'Average edge area =', neMeanRun
write(*,*) 'FPT: treset =', treset
write(*,*)

end program


! initialize RANDOM_SEED
subroutine init_random_seed()
    integer :: values(1:8), k
    integer, dimension(:), allocatable :: seed
    real(8) :: r

    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:) = values(8)
    call random_seed(put=seed)
end subroutine init_random_seed
