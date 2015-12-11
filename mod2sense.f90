module sensing

    use utility

    real(b8), parameter :: g = 1.0
    real(b8), parameter :: g0 = 1000.0
    real(b8), parameter :: gapFlow = 10000.0
    real(b8), parameter :: kappa   = 1.0, mu = 1.0
    real(b8), parameter :: eps = 0.50_b8

contains


! chemical concentration
real(b8) function chemE( x, y)
    implicit none
    real(b8), intent(in) :: x, y
    chemE = g*x + g0
end function chemE


! mean signal in every cell
subroutine getMeanSignal( meanSignal, N, x)
    implicit none
    integer,  intent(in) :: N
    integer,  intent(in),  dimension(:,:,:) :: x
    real(b8), intent(out), dimension(:)     :: meanSignal
    integer :: i, j, nl
    real(b8) :: ds, s, rx, ry

    meanSignal = 0.0
    do i = 1, N
        call occupyCount( nl, x(i,:,:) )
        s = 0.0
        do j = 1, nl
            rx  = real(x(i,j,1))
            ry  = real(x(i,j,2))
            ds = chemE( rx, ry)
            s  = s + ds
        enddo

        meanSignal(i) = s / real(nl)
    enddo

end subroutine getMeanSignal


! measured signal in every cell
subroutine getSpeciesS( meanSignal, N, signal)
    implicit none
    integer,  intent(in) :: N
    real(b8), intent(in),  dimension(:) :: meanSignal
    real(b8), intent(out), dimension(:) :: signal
    real(b8) :: ms
    integer :: i

    signal = 0.0
    do i = 1, N
        ms = meanSignal(i)

        if( ms < 100.0 )then
            call poissonrand( ms, signal(i))
        else
            signal(i) = normal(ms,sqrt(ms))
        endif

        if( signal(i) < 0.0 )then
            signal(i) = 0.0
        endif
    enddo
end subroutine getSpeciesS


! calculate population of x species
subroutine getSpeciesX( N, meanSignal, signal, speciesX)
    implicit none
    integer, intent(in) :: N
    real(b8), intent(in),  dimension(:) :: meanSignal, signal
    real(b8), intent(out), dimension(:) :: speciesX
    real(b8) :: m = 0.0, s = 1.0
    real(b8) :: xi2, xi3, etaX, meanX
    integer :: i

    speciesX = 0.0

    do i = 1, N
        meanX = meanSignal(i)*kappa/mu
        xi2 = normal( m, s)
        xi3 = normal( m, s)
        etaX = sqrt(kappa*meanSignal(i))*xi2 - sqrt(mu*meanX)*xi3

        speciesX(i) = signal(i)*kappa/mu + etaX/mu
        if( speciesX(i) < 0.0 )then
            speciesX(i) = 0.0
        endif
    enddo
end subroutine getSpeciesX


! calculate population of y species
subroutine getSpeciesY( etaY, M, N, signal, y)
    implicit none
    integer, intent(in) :: N
    real(b8),    intent(in),  dimension(:,:) :: M
    real(b8),    intent(in),  dimension(:)   :: etaY, signal
    real(b8),    intent(out), dimension(:)   :: y
    real(b8), allocatable :: b(:), a(:,:)
    integer :: i

    allocate( b(N) )
    allocate( a(N,N) )
    b = kappa * signal + etaY
    a = M

    call gauss_1( a, b, y, N)

    do i = 1, N
        if( y(i) < 0.0 )then
            y(i) = 0.0
        endif
    enddo
    deallocate( a )
    deallocate( b )
end subroutine getSpeciesY


! calculates mean y values for all cells
subroutine getMeanY( meanSignal, M, meanY, N)
    implicit none
    integer, intent(in) :: N
    real(b8),    intent(in),  dimension(:,:) :: M
    real(b8),    intent(in),  dimension(:)   :: meanSignal
    real(b8),    intent(out), dimension(:)   :: meanY
    real(b8), allocatable :: b(:), a(:,:)

    allocate( b(N) )
    allocate( a(N,N) )
    b = kappa * meanSignal
    a = M

    call gauss_1( a, b, meanY, N)

    deallocate( a )
    deallocate( b )
end subroutine


! calculates the noise term in y dynamics
subroutine getEtaY( etaY, gNN, meanS, meanY, N)
    implicit none
    integer, intent(in) :: N
    real(b8),    intent(in),  dimension(:,:) :: gNN
    real(b8),    intent(in),  dimension(:)   :: meanS, meanY
    real(b8),    intent(out), dimension(:)   :: etaY
    real(b8) :: m = 0.0, s = 1.0
    real(b8) :: xi4, xi5, chiJ, sum1
    integer :: i, j

    etaY = 0.0

    do i = 1, N
        xi4 = normal( m, s)
        xi5 = normal( m, s)
        etaY(i) = sqrt( kappa * meanS(i) ) * xi4 - sqrt( mu * meanY(i) ) * xi5

        sum1 = 0.0
        do j = 1, N
            chiJ = normal( m, s)
            sum1 = sum1 + chiJ * sqrt(gNN(i,j)) * ( sqrt(meanY(j)) - sqrt(meanY(i)) )
        enddo

        etaY(i) = etaY(i) + sum1
    enddo

end subroutine getEtaY


! make matrix M for degradation and exchange of y
subroutine makeMtrxM( gNN, M, N)
    implicit none
    integer, intent(in) :: N
    real(b8),    intent(in),  dimension(:,:) :: gNN
    real(b8),    intent(out), dimension(:,:) :: M
    integer :: j, k

    M = 0.0
    do j = 1, N
        M(j,j) = mu + sum( gNN(j,1:N) )
        do k = j+1, N
            if( gNN(j,k) /= 0.0 )then
                M(j,k) = -1.0*gNN(j,k)
                M(k,j) = M(j,k)
            endif
        enddo
    enddo
end subroutine makeMtrxM


! make matri gNN of all gamma_j,k values
subroutine makeMtrxGamma( gNN, N, rSim, sigma, x)
    implicit none
    integer, intent(in) :: N
    real(b8),    intent(inout), dimension(:,:)   :: gNN
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

end subroutine makeMtrxGamma


! calculate cell k's contanct lengths with all its neighbors
subroutine getContactL( k, N, nnLk, rSim, sigma, xcell)
    implicit none
    integer, intent(in) :: N
    integer, intent(out), dimension(:)   :: nnLk
    integer, intent(in),  dimension(:,:) :: sigma
    integer, intent(in),  dimension(:,:) :: xcell
    integer, dimension(2) :: nn, rSim
    integer :: i, inn, j, k, nl

    nnLk = 0

    call occupyCount( nl, xcell)

    do i = 1, nl
        do inn = 1, 4
            call nnGet( inn, nn, rSim, xcell(i,1:2))
            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(nn(1),nn(2)) /= k .AND. sigma(nn(1),nn(2)) /= 0 )then
                nnLk( sigma(nn(1),nn(2))) = nnLk( sigma(nn(1),nn(2))) + 1
            endif
        enddo
    enddo

end subroutine getContactL


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


! returns poisson distributed random number
subroutine poissonrand( lambda, k)
    ! lambda is the mean
    ! k is the random number output
    implicit none
    real(b8), intent(in)  :: lambda
    real(b8), intent(out) :: k
    real(b8) :: l, p, u

    l = exp(-1.0*lambda)
    k = 0.0
    p = 1.0

    do while( p > l )
        k = k + 1.0
        call random_number(u)
        p = p * u
    enddo

    k = k - 1.0

end subroutine poissonrand


subroutine gauss_1(a,b,x,n)
    !============================================================
    ! Solutions to a system of linear equations A*x=b
    ! Method: the basic elimination (simple Gauss elimination)
    ! Alex G. November 2009
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! b(n) - vector of the right hand coefficients b
    ! n - number of equations
    ! output ...
    ! x(n) - solutions
    ! comments ...
    ! the original arrays a(n,n) and b(n) will be destroyed
    ! during the calculation
    !===========================================================
    implicit none
    integer n
    real(b8) a(n,n), b(n), x(n)
    real(b8) c
    integer i, j, k
    !step 1: forward elimination
    !    print*,'hi'
    do k=1, n-1
        do i=k+1,n
            c=a(i,k)/a(k,k)
            a(i,k) = 0.0
            b(i)=b(i)- c*b(k)
            do j=k+1,n
                a(i,j) = a(i,j)-c*a(k,j)
            end do
        end do
    end do
    !step 2: back substitution
    x(n) = b(n)/a(n,n)
    do i=n-1,1,-1
        c=0.0
        do j=i+1,n
            c= c + a(i,j)*x(j)
        end do
        x(i) = (b(i)- c)/a(i,i)
    end do
end subroutine gauss_1


end module
