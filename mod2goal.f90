module goal

use utility

real(b8), parameter :: alpha = 1.0, beta = 2.0
real(b8), parameter :: lambdaA = 1.5, lambdaP = 0.01

contains


! calculate the probability of step being executed
real(b8) function probEval( uNew, uOld, w)
    ! uNew = Energy of potential, new configuration
    ! uOld = Energy of original configuration
    implicit none
    real(b8), intent(in) :: uNew, uOld, w
    real(b8) :: du

    ! calculate the change in the goal function
    du = uNew - uOld

    ! calculate probability
    probEval = exp( min(0.0, -1.0*du + w ) )

end function probEval


! evaluate energy of the configuration
real(b8) function goalEval1( A0, P0, N, rSim, sigma, x)
    ! A0/P0 = preferred area/perimeter of one cell
    ! L = number of lattice sites along one dimension
    ! N = total number of cells
    ! sigma = array of cell labels
    ! x = array of all the cells lattice sites
    implicit none
    real(b8), intent(in) :: P0
    integer,  intent(in) :: A0, N
    integer,  intent(in), dimension(2)     :: rSim
    integer,  intent(in), dimension(:,:)   :: sigma
    integer,  intent(in), dimension(:,:,:) :: x
    integer, dimension(2) :: ls1, ls2
    integer :: i, nx, ny
    real(b8) :: dA, dP, J, sum1, sum2, sum3

    ! sum energy contribution due to J
    sum1 = 0.0
    do nx = 1, rSim(1)
        do ny = 1, rSim(2)

            ls1 = [ nx, ny]

            do i = 1, 2

                call nnGet( i, ls2, rSim, ls1)

                J = jCheck( ls1, ls2, sigma)
                sum1 = sum1 + J

            enddo
        enddo
    enddo

    ! sum energy contribution due to area and perimeter
    sum2 = 0.0
    sum3 = 0.0
    do i = 1, N
        call acCheck( dA, dP, A0, P0, rSim, sigma, x(i,:,:))

        sum2 = sum2 + dA**2.0
        sum3 = sum3 + dP**2.0
    enddo
    sum2 = lambdaA * sum2
    sum3 = lambdaP * sum3

    ! write(*,*) sum1, sum2, sum3
    goalEval1 = sum1 + sum2 + sum3

end function goalEval1


! energy contribution due to J
real(b8) function jCheck( ls1, ls2, sigma)
    ! ls1 = 1st lattice site of interest - x
    ! ls2 = 2nd lattice site of interest - x'
    ! sigma = array of cell labels
    ! alpha = cell-cell energy cost
    ! beta = free-cell energy cost
    implicit none
    integer, dimension(:), intent(in) :: ls1, ls2
    integer, dimension(:,:), intent(in) :: sigma
    integer :: i, j

    i = sigma(ls1(1),ls1(2))

    if( ls2(1) == 0 )then
        j = 0
    else
        j = sigma(ls2(1),ls2(2))
    endif

    if( i == j )then
        jCheck = 0.0
    else
        if( i*j > 0 )then
            jCheck = alpha
        elseif( i*j == 0 )then
            jCheck = beta
        endif
    endif

end function jCheck


! difference in area from preferred area of one cell
real(b8) function aCheck( A0, xcell)
    ! A0 = preferred area of one cell
    ! xcell = array of the lattice points occupied by one cell
    ! a = area of one cell
    implicit none
    integer, intent(in) :: A0
    integer, intent(in), dimension(:,:) :: xcell
    integer :: nl
    real(b8) :: a

    call occupyCount( nl, xcell )
    a = real(nl)

    aCheck = a - real(A0)

end function aCheck


! difference in area from preferred area of one cell
! difference in perimeter from preferred perimeter of one cell
subroutine acCheck( dA, dP, A0, P0, rSim, sigma, xcell)
    ! A0, P0 = preferred area/perimeter of one cell
    ! xcell = array of the lattice points occupied by one cell
    implicit none
    real(b8), intent(in) :: P0
    integer,  intent(in) :: A0
    integer,  intent(in), dimension(:,:) :: sigma
    integer,  intent(in), dimension(:,:) :: xcell
    real(b8), intent(out) :: dA, dP
    integer, dimension(2) :: nn, rSim
    integer :: i, inn, k, nl
    real(b8) :: A, P

    call occupyCount( nl, xcell )
    A = real(nl)

    dA = A - real(A0) ! difference in area

    k = sigma( xcell(1,1), xcell(1,2)) ! cell label

    P = 0.0
    do i = 1, nl
        do inn = 1, 4
            call nnGet( inn, nn, rSim, xcell(i,1:2))
            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(nn(1),nn(2)) /= k )then
                P = P + 1.0
            endif
        enddo
    enddo

    dP = P - P0 ! difference in periemeter

end subroutine acCheck


end module
