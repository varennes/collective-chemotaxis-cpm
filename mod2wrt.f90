module wrtout

use utility

contains


! outputs sigma to fort.101
subroutine wrtSigma( rSim, sigma, tstep)
    ! rlxArea = preferred area of one cell
    ! L = number of lattice sites along one dimension
    ! N = total number of cells
    ! sigma = array of cell labels
    ! x = array of all the cells lattice sites
    implicit none
    integer, intent(in) :: tstep
    integer, intent(in), dimension(2)   :: rSim
    integer, intent(in), dimension(:,:) :: sigma
    integer :: i
    ! write out sigma
    write(101,*) 't =', tstep - 1
    do i = 1, rSim(2) + 2
        write(101,*) sigma(:,i)
    enddo
    write(101,*)

end subroutine wrtSigma


! outputs a, sigma(a), b, sigma(b) to fort.102
subroutine wrtAB( a, b, sigma, tstep)
    implicit none
    integer, intent(in) :: tstep
    integer, dimension(:), intent(in) :: a, b
    integer, dimension(:,:), intent(in) :: sigma

    write(102,*) a, sigma(a(1),a(2)), b, sigma(b(1),b(2)), tstep-1
end subroutine wrtAB


! outputs energy and p, r to fort.104
subroutine wrtU( uNew, uOld, p, r, tstep)
    implicit none
    integer, intent(in) :: tstep
    real(b8), intent(in) ::  uNew, uOld, p, r

    write(104,*) uNew, uOld, p, r, tstep-1
end subroutine wrtU


! outputs x to fort.105
subroutine wrtX( N, x, tstep)
    ! rlxArea = preferred area of one cell
    ! L = number of lattice sites along one dimension
    ! N = total number of cells
    ! sigma = array of cell labels
    ! x = array of all the cells lattice sites
    implicit none
    integer, intent(in) :: N, tstep
    integer, intent(in), dimension(:,:,:) :: x
    integer :: i, j

    ! write out x array
    do i = 1, N
        j = 1
        do while( x(i,j,1) /= 0 )
            write(105,*) x(i,j,1), x(i,j,2), i, tstep-1
            j = j + 1
        enddo
    enddo

end subroutine wrtX


! outputs edge to fort.103
subroutine wrtEdge( edge, rSim, sigma, tstep)
    ! L = number of lattice sites along one dimension
    ! sigma = array of cell labels
    implicit none
    integer, intent(in) :: tstep
    integer, intent(in), dimension(2)   :: rSim
    integer, intent(in), dimension(:,:) :: sigma
    integer, intent(in), dimension(:,:) :: edge
    integer, dimension(rSim(1)+2,rSim(2)+2) :: sigmaWrt
    integer :: i, j, k, n

    sigmaWrt = sigma

    n = 1
    do while( edge(n,1) /= 0 )
        i = edge(n,1)
        j = edge(n,2)
        sigmaWrt(i,j) = 9
        n = n + 1
    enddo

    write(103,*) 't =', tstep - 1
    do i = 1, rSim(2) + 2
        write(103,*) sigmaWrt(:,i)
    enddo
    write(103,*)

end subroutine wrtEdge


subroutine wrtEdgeArray( edge, tstep)
    implicit none
    integer, intent(in) :: tstep
    integer, intent(in), dimension(:,:) :: edge
    integer :: i

    i = 1
    do while( edge(i,1) /= 0 )
        write(106,*) edge(i,:), tstep
        i = i + 1
    enddo

end subroutine wrtEdgeArray


subroutine wrtPolar( N, p, tstep)
    implicit none
    integer, intent(in) :: N, tstep
    real(b8),    intent(in), dimension(:,:) :: p
    integer :: i

    do i = 1, N
        write(140,*) p(i,:), tstep - 1
    enddo

end subroutine wrtPolar

end module
