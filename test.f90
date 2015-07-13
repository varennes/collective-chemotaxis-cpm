program test1

use utility
use polar
use sensing
use wrtout

implicit none

integer :: i, j, ne, ne0, nf, nl, check, aSig, bSig
integer :: A0, N, x1, x2
integer, dimension(2) :: r0, rCell, rSim

integer, allocatable :: sigma(:,:), x(:,:,:), nnL(:,:)

real(b8) :: rx1, rx2
real(b8), allocatable :: cellCOM(:,:), speciesR(:), speciesX(:), speciesY(:)
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

allocate( nnL(N,N) )
allocate(    sigma( rSim(1) + 2, rSim(2) + 2) )
allocate(     x( N, 4*A0, 2) )
allocate( meanSignal(N) )
allocate( speciesR(N) )
allocate( speciesX(N) )
allocate( speciesY(N) )
allocate( meanY(N) )
allocate(  etaY(N) )
allocate( signal(N) )
allocate( gNN(N,N) )
allocate(   M(N,N) )

allocate( cellCOM(N,2) )

x    = 0
sigma    = 0

! initialize sigma
rSim(1) = x1
call itlSigma( r0, rCell, rSim, sigma)

! call wrtSigma( rSim, sigma, 1)

! initialize x array
rSim(1) = x1
call makeX( N, rSim, sigma, x)

! calculate contact lengths
do i = 1, N
    call getContactL( i, N, nnL(i,:), rSim, sigma, x(i,:,:))
enddo

do i = 1, N
    call calcCellCOM( x(i,:,:),  cellCOM(i,:))
    write(155,*) cellCOM(i,:), i
enddo

do i = 2, rSim(1)+1
    rx1 = real(i)
    rx2 = real(rSim(2))
    write(122,*) chemE( rx1, rx2)
enddo

call getMeanSignal( meanSignal, N, x)
call getSpeciesS( meanSignal, N, signal)

call getSpeciesX( N, meanSignal, signal, speciesX)

call makeMtrxGamma( gNN, N, rSim, sigma, x)
call makeMtrxM( gNN, M, N)

call getMeanY( meanSignal, M, meanY, N)
call getEtaY( etaY, gNN, meanSignal, meanY, N)

call getSpeciesY( etaY, M, N, signal, speciesY)

! do i = 1, N
!     write(133,*) gNN(i,:)
!     write(144,*)   M(i,:)
! enddo


do i = 1, 1000
    call getSpeciesS( meanSignal, N, signal)
    call getSpeciesX( N, meanSignal, signal, speciesX)
    call getEtaY( etaY, gNN, meanSignal, meanY, N)
    call getSpeciesY( etaY, M, N, signal, speciesY)

    speciesR = speciesX - speciesY

    write(111,*) signal(:)
    write(112,*) speciesX(:)
    write(113,*) speciesY(:)
    write(114,*) speciesR(:)
enddo

write(*,*) 'N =', N
! write(*,*) '      xCOM =', xCOM
write(*,*) 'meanSignal =', meanSignal
write(*,*)
write(*,*) '  contact lengths'
do i = 1, N
    write(*,*) ' cell',i,'nnL =',nnL(i,:)
enddo

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
