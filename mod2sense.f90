module sensing

    use utility

    ! b8 will be used to define reals with 14 digits
    integer, parameter:: b8 = selected_real_kind(14)
    real, parameter :: g = 0.0025000000
    real, parameter :: gapFlow = 1.00000

    real, parameter :: kappa = 1.0000, mu = 1.0000

contains


! chemical concentration
real function chemE( x, y)
    implicit none
    real, intent(in) :: x, y
    chemE = g*x + 1.05000
end function chemE


! mean signal in a cell
real function getMeanSignal( xcell)
    implicit none
    integer, intent(in), dimension(:,:) :: xcell
    integer :: i, nl
    real :: ds, s, x, y

    call occupyCount( nl, xcell )

    s = 0.000000
    do i = 1, nl
        x  = real(xcell(i,1))
        y  = real(xcell(i,2))
        ds = chemE( x, y)
        s  = s + ds
    enddo

    getMeanSignal = s / real(nl)

end function getMeanSignal


! measured signal in a cell
real function getCellSignal( meanSignal)
    implicit none
    real, intent(in) :: meanSignal
    real(b8) :: ms

    ms = meanSignal
    getCellSignal = normal(ms,ms)
    if( getCellSignal < 0.0 )then
        getCellSignal = 0.00000
    endif
end function getCellSignal


! local population of x species
real function getLocalX( meanSignal, signal)
    implicit none
    real, intent(in) :: meanSignal, signal
    real(b8) :: m = 0.0, s = 1.0
    real :: eta2, eta3, etaX, meanX

    meanX = meanSignal*kappa/mu
    eta2 = normal( m, s)
    eta3 = normal( m, s)
    etaX = sqrt(kappa*meanSignal)*eta2 - sqrt(mu*meanX)*eta3

    getLocalX = signal*kappa/mu + etaX/mu
    if( getLocalX < 0.0 )then
        getLocalX = 0.0
    endif
end function getLocalX


subroutine gammaMat( gNN, N, rSim, sigma, x)
    implicit none
    integer, intent(in) :: N
    real,    intent(inout), dimension(:,:)   :: gNN
    integer, intent(in),    dimension(:,:)   :: sigma
    integer, intent(in),    dimension(:,:,:) :: x
    integer, dimension(2) :: nn, rSim
    integer, allocatable  :: neighborL(:)
    integer :: i, inn, j, k, nl

    allocate( neighborL(N) )

    gNN = 0.0
    do j = 1, N

        neighborL = 0
        call occupyCount( nl, x(j,:,:) )

        do i = 1, nl
            do inn = 1, 4
                call nnGet( inn, nn, rSim, x(j,i,1:2))
                if( nn(1) == 0 )then
                    cycle
                endif
                if( sigma(nn(1),nn(2)) /= j .AND. sigma(nn(1),nn(2)) /= 0 )then
                    neighborL( sigma(nn(1),nn(2))) = neighborL( sigma(nn(1),nn(2))) + 1
                endif
            enddo
        enddo

        do k = j+1, N
            if ( neighborL(k) == 0 ) then
                cycle
            end if

            gNN(j,k) = gapFlow * real(neighborL(k))
            gNN(k,j) = gNN(j,k)
        enddo

    enddo

end subroutine gammaMat


! returns random number between 0 - 1
function ran1()
    implicit none
    real(b8) ran1,x
    call random_number(x) ! built in fortran 90 random number function
    ran1=x
end function ran1


! returns a normal distribution
function normal(mean,sigma)
    implicit none
    real(b8) normal,tmp
    real(b8) mean,sigma
    integer flag
    real(b8) fac,gsave,rsq,r1,r2
    save flag,gsave
    data flag /0/
    if (flag.eq.0) then
    rsq=2.0_b8
        do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do
            r1=2.0_b8*ran1()-1.0_b8
            r2=2.0_b8*ran1()-1.0_b8
            rsq=r1*r1+r2*r2
        enddo
        fac=sqrt(-2.0_b8*log(rsq)/rsq)
        gsave=r1*fac
        tmp=r2*fac
        flag=1
    else
        tmp=gsave
        flag=0
    endif
    normal=tmp*sigma+mean
    return
end function normal


end module
