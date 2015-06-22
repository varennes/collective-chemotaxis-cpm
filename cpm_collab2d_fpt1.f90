program cpmcollab

! boundary condition:  If pickb chooses a lattice point outside of the space
!                      pick again without advancing the clock.

use utility
use goal
use simpleconnect
use wrtout

! allocate variables
implicit none
integer :: i, j, ne, ne0, nf, nl, check
integer :: tcount, tELEM, tMCS, tmax
integer :: nRun, runTotal

integer :: A0, N, x1, x2
integer, dimension(2) :: r0, rCell, rSim

integer, allocatable :: sigma(:,:), sigmaTmp(:,:)
integer, allocatable :: x(:,:,:), xTmp(:,:,:)
integer, allocatable :: edge(:,:), filled(:,:), node(:)
integer, dimension(2) :: a, b, nn

real, allocatable :: firstpass(:), xCOM(:,:)
real :: d, df, p, r, uNew, uOld
real :: neMean, neMeanRun
real :: t0, tf

call cpu_time(t0)

! initialize variables
open(unit=11,file='input.txt',status='old',action='read')
read(11,*) r0(1), r0(2)
read(11,*) rCell(1), rCell(2)
read(11,*) x1, x2
read(11,*) rSim(2)
read(11,*) tmax
read(11,*) runTotal
read(11,*) df

rSim(1) = x1 + x2
N       = rCell(1) * rCell(2)
A0      = r0(1) * r0(2)

write(*,*) '   N =',N
write(*,*) '  A0 =',A0
write(*,*) 'rSim =',rSim
write(*,*) '  df =',df
write(*,*)

call init_random_seed()

allocate(    sigma( rSim(1), rSim(2)) )
allocate( sigmaTmp( rSim(1), rSim(2)) )
allocate(     x( N, 2*A0, 2) )
allocate(  xTmp( N, 2*A0, 2) )
allocate(  xCOM( tmax, 2) )

allocate( edge( rSim(1)*rSim(2), 2) )
allocate( filled( 2*A0, 2) )
allocate( node(2) )

allocate( firstpass(runTotal) )

neMeanrun = 0.0
firstpass = 0.0

do nRun = 1, runTotal

tcount = 1
tELEM  = 1 ! elementary time step
tMCS   = 1 ! MCS step, 1 MCS step = L^2 ELEM steps

edge   = 0
node   = 0
filled = 0

x    = 0
xTmp = 0
sigma    = 0
sigmaTmp = 0

d    = 0.0
xCOM = 0.0

a = 1
b = 1

! initialize sigma
rSim(1) = x1
call itlSigma( r0, rCell, rSim, sigma)

! initialize edge
call itlEdge( edge, ne, rSim, sigma)
ne0 = ne
neMean = real(ne)

! initialize x array
rSim(1) = x1
call makeX( N, rSim, sigma, x)

! initialize xCOM
call calcXCOM( N, x, xCOM(1,:))

! calculate initial energy
rSim(1) = x1 + x2
uOld = goalEval1( A0, N, rSim, sigma, x)
uNew = 0.0

! write outputs
! write(*,*) tELEM,' xCOM =',xCOM(tELEM,:)
! call wrtSigma( rSim, sigma, tELEM)
! call wrtEdge( edge, rSim, sigma, tELEM)
! call wrtEdgeArray( edge, tELEM)
! call wrtU( 0.0, uOld, 0.0, 0.0, tELEM)
! call wrtX( N, x, tELEM)

do while( tMCS < tmax )
    tELEM = tELEM + 1

    call countEdge( edge, ne)

    call picka( a, edge, ne)

    call pickb( a, b, rSim)

    if( sigma(b(1),b(2)) /= sigma(a(1),a(2)) )then
        ! make and update sigmaTmp
        sigmaTmp = sigma
        ! sigmaTmp(b(1),b(2)) = aSigma
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
            uNew = goalEval1( A0, N, rSim, sigmaTmp, xTmp)
            p    = probEval( uNew, uOld)

            call random_number(r)

            if( r < p )then
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

        ! calculate COM
        call calcXCOM( N, x, xCOM(tMCS,:))

        ! write outputs
        ! write(*,*) tMCS,' xCOM =',xCOM(tMCS,:)
        ! call wrtSigma( rSim, sigma, tMCS)
        ! call wrtX( N, x, tMCS)

        ! calculate d
        d = calcD( xCOM(tMCS,1), xCOM(1,1))
        if( d >= df )then
            firstpass(nRun) = tMCS - 1
            tMCS = tmax
        elseif( d <= -1.0*df )then
            firstpass(nRun) = 0.0
            tMCS = tmax
        endif

    endif

enddo ! end while loop

neMean    = neMean / real(tcount)
neMeanRun = neMeanRun + neMean

enddo ! end run loop

do i = 1, runTotal
    write(109,*) firstpass(i), i
enddo

neMeanRun = neMeanRun / real(runTotal)

close(11)

call cpu_time(tf)
write(*,*) 'Run time =',tf-t0
write(*,*) 'Initial edge area =',ne0
write(*,*) 'Average edge area =', neMeanRun

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
