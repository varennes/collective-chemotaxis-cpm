module goal

use utility

contains


! calculate the probability of step being executed
real function probEval( uNew, uOld, w)
    ! uNew = Energy of potential, new configuration
    ! uOld = Energy of original configuration
    implicit none
    real, intent(in) :: uNew, uOld, w
    real :: du

    ! calculate the change in the goal function
    du = uNew - uOld

    ! calculate probability
    probEval = exp( min(0.0, -1.0*du + w ) )

end function probEval


! evaluate energy of the configuration
real function goalEval1( A0, N, rSim, sigma, x)
    ! A0 = preferred area of one cell
    ! L = number of lattice sites along one dimension
    ! N = total number of cells
    ! sigma = array of cell labels
    ! x = array of all the cells lattice sites
    implicit none
    integer, intent(in) :: A0, N
    integer, intent(in), dimension(2)     :: rSim
    integer, intent(in), dimension(:,:)   :: sigma
    integer, intent(in), dimension(:,:,:) :: x
    integer, dimension(2) :: ls1, ls2
    integer :: i, nx, ny
    real :: dA, J, lambda, sum1, sum2

    lambda = 1.0

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

    ! sum energy contribution due to area
    sum2 = 0.0
    do i = 1, N
        dA = aCheck( A0, x(i,:,:))

        sum2 = sum2 + dA**2.0
    enddo
    sum2 = lambda * sum2

    goalEval1 = sum1 + sum2

end function goalEval1


! energy contribution due to J
real function jCheck( ls1, ls2, sigma)
    ! ls1 = 1st lattice site of interest - x
    ! ls2 = 2nd lattice site of interest - x'
    ! sigma = array of cell labels
    ! alpha = cell-cell energy cost
    ! beta = free-cell energy cost
    implicit none
    integer, dimension(:), intent(in) :: ls1, ls2
    integer, dimension(:,:), intent(in) :: sigma
    integer :: i, j
    real :: alpha, beta

    alpha = 0.5
    beta  = 1.0

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
real function aCheck( A0, xcell)
    ! A0 = preferred area of one cell
    ! xcell = array of the lattice points occupied by one cell
    ! a = area of one cell
    implicit none
    integer, intent(in) :: A0
    integer, intent(in), dimension(:,:) :: xcell
    integer :: nl
    real :: a

    call occupyCount( nl, xcell )
    a = real(nl)

    aCheck = a - real(A0)

end function aCheck


end module
