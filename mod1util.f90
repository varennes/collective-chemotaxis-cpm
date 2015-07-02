module utility

contains


! create an array of all the cells lattice sites
subroutine makeX( N, rSim, sigma, x)
    ! L = number of lattice sites along one dimension
    ! N = total number of cells
    ! sigma = array of cell labels
    ! x = array of all the cells lattice sites
    ! lsCount = array that keeps track the number of lattice sites
    implicit none
    integer, intent(in) :: N
    integer, intent(in),  dimension(2)     :: rSim
    integer, intent(in),  dimension(:,:)   :: sigma
    integer, intent(out), dimension(:,:,:) :: x
    integer :: cellIndex, i, nx, ny
    integer, allocatable :: lsCount(:)

    allocate( lsCount(N) )
    lsCount(:) = 0
    x = 0

    do nx = 1, rSim(1) + 2
        do ny = 1, rSim(2) + 2
            if( sigma(nx,ny) /= 0 )then

                cellIndex = sigma(nx,ny)
                lsCount(cellIndex) = lsCount(cellIndex) + 1

                i = lsCount(cellIndex)
                x(cellIndex,i,:) = [nx,ny]
            endif
        enddo
    enddo

    deallocate( lsCount )
end subroutine makeX


! initialize the array sigma
subroutine itlSigma( r0, rCell, rSim, sigma)
    ! r0 = initial dimnesions of individual cells
    ! rCell = cluster dimensions in terms of number of cells
    ! rSim = dimensions of cell initialization space
    ! sigma = array of cell labels
    implicit none
    integer, intent(in),    dimension(2)   :: r0, rCell, rSim
    integer, intent(inout), dimension(:,:) :: sigma
    integer, dimension(2) :: rGroup
    integer :: i, nx, ny, xi, xf, yi, yf

    sigma = 0

    rGroup = [ r0(1)*rCell(1), r0(2)*rCell(2) ]

    if( rGroup(1) > rSim(1) .OR. rGroup(2) > rSim(2) )then
        write(*,*) 'WARNING: Not enough of room for cell cluster!'
    endif

    xi = ceiling( real(rSim(1))/2.0 - real(rGroup(1))/2.0 ) + 2

    i = 0
    do nx = 1, rCell(1)

        xf = xi + r0(1) - 1
        yi = ceiling( real(rSim(2))/2.0 - real(rGroup(2))/2.0 ) + 2

        do ny = 1, rCell(2)

            yf = yi + r0(2) - 1

            i = i + 1
            sigma( xi:xf, yi:yf) = i

            yi = yf + 1
        enddo

        xi = xf + 1
    enddo

end subroutine itlSigma


subroutine itlEdge( edge, ne, rSim, sigma)
    ! edge = array of lattice sites with diff. sigma neighbors
    implicit none
    integer, intent(out) :: ne
    integer, intent(out), dimension(:,:) :: edge
    integer, intent(in),  dimension(2)   :: rSim
    integer, intent(in),  dimension(:,:) :: sigma
    integer, dimension(2) :: nn
    integer :: i1, i2, j, k

    ne = 1
    edge = 0

    do i1 = 1, rSim(1) + 2
    do i2 = 1, rSim(2) + 2
        do k = 1, 4
            call nnGet( k, nn, rSim, [i1,i2])

            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(nn(1),nn(2)) /= sigma(i1,i2) )then
                edge(ne,:) = [i1,i2]
                ne = ne + 1
                exit
            endif
        enddo
    enddo
    enddo
    ne = ne - 1

end subroutine itlEdge


! calculate center of mass of the group of cells
subroutine calcXCOM( N, x, xCOMt)
    ! N = total number of cells
    ! x = array of all the cells lattice sites
    ! xCOMt = array of center of mass of all cells at some time
    implicit none
    integer, intent(in) :: N
    integer, intent(in),  dimension(:,:,:) :: x
    real,    intent(out), dimension(:)     :: xCOMt
    integer :: i, j, k, nl, nlSum

    nlSum = 0
    xCOMt = 0.0

    do i = 1, N
        call occupyCount( nl, x(i,:,:))

        do j = 1, nl
            do k = 1, 2
                xCOMt(k) = xCOMt(k) + real(x(i,j,k))
            enddo
        enddo

        nlSum = nlSum + nl
    enddo

    xCOMt = xCOMt / real(nlSum)

end subroutine calcXCOM


! calculate displacement along one dimension
real function calcD( xCOMt, x0)
    ! xCOMt = COM position of cell group at a pre-set time
    ! x0 = initial COM position of cell group
    implicit none
    real, intent(in) :: xCOMt, x0

    calcD = xCOMt - x0

end function calcD


! calculate displacement along one dimension
real function calcMSD( xCOMt, x0)
    ! xCOMt = COM array of cell group at a pre-set time
    ! x0 = initial COM array of cell group
    implicit none
    real, intent(in), dimension(2) :: xCOMt, x0

    calcMSD = (xCOMt(1)-x0(1))**2.0 + (xCOMt(2)-x0(2))**2.0

end function calcMSD


! Output coordinates of the nearest neighbor (nn)
subroutine nnGet( i, nn, rSim, x)
    ! i indicates the nn we are interested in
    ! i = 1 -- nn up
    ! i = 2 -- nn right
    ! i = 3 -- nn down
    ! i = 4 -- nn left
    ! nn = array of nearest neighbor coordinates
    ! rSim = dimensions of simulation space
    ! x = coordinates of point
    implicit none
    integer, intent(in) :: i
    integer, dimension(2), intent(in) :: rSim, x
    integer, dimension(2), intent(out) :: nn

    if( i == 1 )then
        nn(1) = x(1) + 1
        nn(2) = x(2)

    elseif( i == 2 )then
        nn(1) = x(1)
        nn(2) = x(2) + 1

    elseif( i == 3 )then
        nn(1) = x(1) - 1
        nn(2) = x(2)

    elseif( i == 4 )then
        nn(1) = x(1)
        nn(2) = x(2) - 1

    endif

    if( nn(1) > (rSim(1) + 2) .OR. nn(1) < 1 )then
        nn = [ 0, 0]
    elseif( nn(2) > (rSim(2) + 2) .OR. nn(2) < 1 )then
        nn = [ 0, 0]
    endif

end subroutine nnGet


! Count the number of occupied lattice points
subroutine occupyCount( nl, xcell )
    ! xcell =  array of lattice sites occupied by one cell x(i,:,:)
    ! nl = number of lattice sites contained in xcell
    implicit none
    integer, dimension(:,:), intent(in) :: xcell
    integer, intent(out) :: nl

    ! count how many lattice sites one cell contains
    nl = 1
    do while( xcell(nl,1) /= 0 )
        nl = nl + 1
    enddo
    nl = nl - 1
end subroutine occupyCount


end module
