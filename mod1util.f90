module utility

    ! b8 will be used to define reals with 14 digits
    integer, parameter:: b8 = selected_real_kind(14)

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


subroutine itlSigmaRandom( N, r0, rCell, rSim, sigma)
    ! r0 = initial dimnesions of individual cells
    ! rCell = cluster dimensions in terms of number of cells
    ! rSim = dimensions of cell initialization space
    ! sigma = array of cell labels
    implicit none
    integer, intent(in) :: N
    integer, intent(in),    dimension(2)   :: r0, rCell, rSim
    integer, intent(inout), dimension(:,:) :: sigma
    integer, allocatable :: cellGrid(:,:), ocpyGrid(:,:)
    integer, dimension(4,2) :: dr
    integer, dimension(2) :: center, cell1, cellk
    integer :: count, i, j, k
    real(b8) :: r

    allocate( cellGrid(N,2) )
    allocate( ocpyGrid(N+1,N+1) )

    sigma    = 0
    cellGrid = 0
    ocpyGrid = 0
    ocpyGrid(N/2+1,N/2+1) = 1
    cellGrid(1,:) = [ N/2+1, N/2+1]

    dr(1,:) = [ 1, 0]
    dr(2,:) = [ 0, 1]
    dr(3,:) = [-1, 0]
    dr(4,:) = [ 0,-1]

    center(1) = rSim(1)/2 + 2
    center(2) = rSim(2)/2 + 2
    cell1(1)  = center(1) - r0(1)/2
    cell1(2)  = center(2) - r0(2)/2

    sigma( cell1(1):cell1(1)+r0(1)-1, cell1(2):cell1(2)+r0(2)-1) = 1

    i = 1
    k = 1
    do i = 1, N-1
        if( k >= N )then
            exit
        endif

        call random_number(r)
        j = 1 + floor(r*4)

        do count = 1, 4
            if( k >= N )then
                exit
            endif

            if( ocpyGrid( cellGrid(i,1)+dr(j,1), cellGrid(i,2)+dr(j,2)) == 0 )then
                k = k + 1
                cellGrid(k,1) = cellGrid(i,1)+dr(j,1)
                cellGrid(k,2) = cellGrid(i,2)+dr(j,2)

                ocpyGrid( cellGrid(k,1), cellGrid(k,2)) = k
            endif

            j = j + 1
            if( j > 4 )then
                j = 1
            endif
        enddo
    enddo

    do i = 1, N
        write(177,*) ocpyGrid(i,:)
    enddo

    do k = 2, N
        i = cellGrid(1,1) - cellGrid(k,1) ! grid space distance from cell 1
        j = cellGrid(1,1) - cellGrid(k,2)

        cellk(1) = cell1(1) + i*r0(1)
        cellk(2) = cell1(2) + j*r0(2)

        if( cellk(1) < 2 .OR. (cellk(1)+r0(1)-1) > rSim(1)+1 )then
            write(*,*) 'Cells out of bounds! You done goofed!'
            exit
        endif
        if( cellk(2) < 2 .OR. (cellk(2)+r0(1)-1) > rSim(2)+1 )then
            write(*,*) 'Cells out of bounds! You done goofed!'
            exit
        endif

        sigma( cellk(1):cellk(1)+r0(1)-1, cellk(2):cellk(2)+r0(2)-1) = k
    enddo

end subroutine itlSigmaRandom


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
    real(b8),    intent(out), dimension(:)     :: xCOMt
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
real(b8) function calcD( xCOMt, x0)
    ! xCOMt = COM position of cell group at a pre-set time
    ! x0 = initial COM position of cell group
    implicit none
    real(b8), intent(in) :: xCOMt, x0

    calcD = xCOMt - x0

end function calcD


! calculate displacement along one dimension
real(b8) function calcMSD( xCOMt, x0)
    ! xCOMt = COM array of cell group at a pre-set time
    ! x0 = initial COM array of cell group
    implicit none
    real(b8), intent(in), dimension(2) :: xCOMt, x0

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
