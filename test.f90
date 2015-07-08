program test1

use utility
use sensing
use wrtout

implicit none

integer :: i, j, ne, ne0, nf, nl, check, aSig, bSig
integer :: A0, N, x1, x2
integer, dimension(2) :: r0, rCell, rSim

integer, allocatable :: sigma(:,:), x(:,:,:)

real(b8) :: rx1, rx2
real(b8), allocatable :: y(:)
real(b8), allocatable :: etaY(:), gNN(:,:), M(:,:), meanSignal(:), signal(:), meanY(:)
real(b8), dimension(2) :: xCOM

! initialize parameters, variables
open(unit=11,file='input.txt',status='old',action='read')
read(11,*) r0(1), r0(2)
read(11,*) rCell(1), rCell(2)
read(11,*) x1, x2
read(11,*) rSim(2)
! read(11,*) tmax
! read(11,*) runTotal
! read(11,*) df

rSim(1) = x1 + x2
N       = rCell(1) * rCell(2)
A0      = r0(1) * r0(2)

! initialize parameters for polarization
! open(unit=12,file='polarInput.txt',status='old',action='read')
! read(12,*) plrP
! read(12,*) plrR

write(*,*) '   N =',N
write(*,*) '  A0 =',A0
write(*,*) 'rSim =',rSim
! write(*,*) '  df =',df
! write(*,*) 'plrP =',plrP,' plrR =',plrR
write(*,*)

call init_random_seed()

allocate(    sigma( rSim(1) + 2, rSim(2) + 2) )
allocate(     x( N, 4*A0, 2) )
allocate( meanSignal(N) )
allocate( y(N) )
allocate( meanY(N) )
allocate(  etaY(N) )
allocate( signal(N) )
allocate( gNN(N,N) )
allocate(   M(N,N) )

x    = 0
sigma    = 0

! initialize sigma
rSim(1) = x1
call itlSigma( r0, rCell, rSim, sigma)

call wrtSigma( rSim, sigma, 1)

! initialize x array
rSim(1) = x1
call makeX( N, rSim, sigma, x)

do i = 1, N
    meanSignal(i) = getMeanSignal( x(i,:,:))
    signal(i) = getCellSignal( meanSignal(i))
enddo


! do i = 1, 1000
!     write(111,*) getCellSignal( meanSignal(1))
!     write(112,*) getLocalX( meanSignal(1), signal(1))
! enddo

do i = 2, rSim(1)+1
    rx1 = real(i)
    rx2 = real(rSim(2))
    write(122,*) chemE( rx1, rx2)
enddo

call makeMtrxGamma( gNN, N, rSim, sigma, x)
call makeMtrxM( gNN, M, N)
call getMeanY( meanSignal, M, meanY, N)

do i = 1, N
    write(133,*) gNN(i,:)
    write(144,*)   M(i,:)
enddo

call calcXCOM( N, x, xCOM)

call getEtaY( etaY, gNN, meanSignal, meanY, N)

call getSpeciesY( etaY, M, N, signal, y)

do i = 1, 1000
    write(111,*) getCellSignal( meanSignal(1))
    write(112,*) getLocalX( meanSignal(1), signal(1))
    call getEtaY( etaY, gNN, meanSignal, meanY, N)
    call getSpeciesY( etaY, M, N, signal, y)
    write(113,*) y(1)
enddo

write(*,*) 'N =', N
write(*,*) '      xCOM =', xCOM
write(*,*) 'meanSignal =', meanSignal
write(*,*)
write(*,*) '     meanY =', meanY
write(*,*) '         y =', y

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
