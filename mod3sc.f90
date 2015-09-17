module simpleconnect

use utility

contains


! count number of non-zero values in array edge
subroutine countEdge( edge, ne)
    implicit none
    integer, intent(out) :: ne
    integer, dimension(:,:), intent(in) :: edge

    ne = 1
    do while( edge(ne,1) /= 0 )
        ne = ne + 1
    enddo
    ne = ne - 1
end subroutine countEdge


! pick a from edge
subroutine picka( a, edge, ne)
    implicit none
    integer, intent(in) :: ne
    integer, intent(in),  dimension(:,:) :: edge
    integer, intent(out), dimension(2)   :: a
    integer :: na
    real(b8) :: r

    call random_number(r)
    na = 1 + floor( real(ne)*r )

    ! find out which lattice point edge(na) corresponds to
    a = edge(na,:)

end subroutine picka


! pick b
subroutine pickb( a, b, rSim)
    implicit none
    integer, intent(in),  dimension(2) :: a, rSim
    integer, intent(out), dimension(2) :: b
    integer :: i, j
    real(b8) :: r

    b = 0
    do while( b(1) == 0 )
        call random_number(r)
        i = 1 + floor(4.0*r)
        call nnGet( i, b, rSim, a)
    enddo

end subroutine pickb


! update array edge. Check whether lattice site ls is an edge.
subroutine updateEdge( ls, edge, ne, rSim, sigma)
    implicit none
    integer, intent(inout) :: ne
    integer, dimension(2), intent(in) :: ls, rSim
    integer, dimension(:,:), intent(inout) :: edge
    integer, dimension(:,:), intent(in) :: sigma
    integer, dimension(2) :: nn
    integer :: check, i, j

    check = 0

    ! check if ls is an edge
    do i = 1, 4

        call nnGet( i, nn, rSim, ls)
        if( nn(1) == 0 )then
            cycle
        endif

        if( sigma(ls(1),ls(2)) /= sigma(nn(1),nn(2)) )then
            ! ls is an edge
            check = 1

            ! check whether ls is in the array edge
            do j = 1, ne
                if( edge(j,1) == ls(1) .AND. edge(j,2) == ls(2) )then
                    ! ls is in array edge
                    exit
                elseif( j == ne )then
                    if( edge(j,1) /= ls(1) .AND. edge(j,2) /= ls(2) )then
                        ! ls is not in the array edge
                        ne = ne + 1
                        edge(ne,:) = ls
                        exit
                    endif
                endif
            enddo

            exit
        endif
    enddo

    if( check == 0 )then
        ! ls is not an edge
        do i = 1, ne
            if( edge(i,1) == ls(1) .AND. edge(i,2) == ls(2) )then
                call deleteEdge( i, edge, ne)
            endif
        enddo
    endif

end subroutine updateEdge


subroutine deleteEdge( i, edge, ne)
    implicit none
    integer, intent(in) :: i
    integer, intent(inout) :: ne
    integer, dimension(:,:), intent(inout) :: edge

    if( i == ne )then
        edge(i,:) = [ 0, 0]
    else
        edge( i : ne - 1, :) = edge( i + 1 : ne, : )
        edge(ne,:) = [ 0, 0]
    endif

    ne = ne - 1
end subroutine deleteEdge


! Checks whether a cell is simply connected or not using flood fill algorithm
recursive subroutine floodFill( node, filled, rSim, xcell)
    ! L = number of lattice sites along one dimension
    ! node = array of lattice site coordinates of node
    ! filled = array of all lattice sites that have been filled by the algorithm
    ! xcell =  array of lattice sites occupied by cell x(i,:,:)
    implicit none
    integer, dimension(2), intent(in) :: node, rSim
    integer, dimension(:,:), intent(inout) :: filled
    integer, dimension(:,:), intent(in) :: xcell
    integer, dimension(2) :: nn
    integer :: i, j, nf, nl

    call occupyCount( nf, filled )
    call occupyCount( nl, xcell)

    do i = 1, nf
        if( node(1) == filled(i,1) .AND. node(2) == filled(i,2) )then
            return
        endif
    enddo

    j = 0
    do i = 1, nl
        if( node(1) == xcell(i,1) .AND. node(2) == xcell(i,2) )then
            j = 1
        endif
    enddo
    if( j /= 1 )then
        return
    endif

    nf = nf + 1
    if( nf > nl )then
        return
    endif

    filled( nf, :) = node

    call nnGet( 1, nn, rSim, node)
    call floodFill( nn, filled, rSim, xcell)

    call nnGet( 2, nn, rSim, node)
    call floodFill( nn, filled, rSim, xcell)

    call nnGet( 3, nn, rSim, node)
    call floodFill( nn, filled, rSim, xcell)

    call nnGet( 4, nn, rSim, node)
    call floodFill( nn, filled, rSim, xcell)

end subroutine floodFill


end module
