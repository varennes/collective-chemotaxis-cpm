module wrtout

use utility

contains


! calculate the mean cluster size in the system
! MUST HAVE UP TO DATE nnL array. run "getContactL" to update nnL.
! output mean cluster size to fort.230+nRun
! output list of different clusters at a given timestep to 210+nRun
subroutine wrtClstrSize( N, nnL, nRun, tMCS)
    implicit none
    integer, intent(in) :: N, nRun, tMCS
    integer, intent(in), dimension(:,:) :: nnL
    integer, dimension(N+2,N+2) :: clusterList
    integer :: i, ii, in, j, k, l, ll
    integer :: ic, icTotal, clCheck, nc, ncTotal, ncNew, ncOld
    real(b8) :: nsub
    ! output cluster sizes and mean subcluster size
    ic = 1
    icTotal = 1
    ncTotal = 1
    clusterList(:,:) = 0
    clusterList(icTotal,1) = 1
    do i = 1, N
        clCheck = 0
        do ii = 1, icTotal
            ! do in = 1, ncTotal
            do in = 1, N
                if( clusterList(ii,in) == i )then
                    clCheck = 1
                    exit
                endif
            enddo
            if( clCheck == 1 )then
                exit
            endif
        enddo
        if( clCheck == 0 )then
            icTotal = icTotal + 1
            clusterList(icTotal,1) = i
            ic = icTotal
        endif
        ncTotal = 0
        do in = 1, N
            if( clusterList(ic,in) /= 0 )then
                ncTotal = ncTotal + 1
            endif
        enddo
        ncOld = 0
        ncNew = ncTotal
        ! write(*,*) ''
        ! write(*,*) 'i = ', i,'| icTotal = ',icTotal,' ncTotal= ',ncTotal
        do while( ncOld /= ncNew )
            nc = 0
            ncOld = ncNew
            ! write(*,*)
            ! write(*,*) 'while strt, ncOld =',ncOld
            do j = 1, ncOld
                k = clusterList(ic,j)
                do l = 1, N
                    clCheck = 0
                    if( nnL(k,l) /= 0 )then
                        ! write(*,*) 'nnL /= 0 for k=',k, '& l=',l,' clCheck=',clCheck
                        do ll = 1, ncOld+nc
                            if( clusterList(ic,ll) == l )then
                                clCheck = 1
                                exit
                            endif
                        enddo
                        if( clCheck == 0 )then
                            nc = nc + 1
                            clusterList(ic,ncOld+nc) = l
                            ! write(*,*) 'ncOld+nc =',ncOld+nc,'l=',l
                        endif
                    endif
                enddo
            enddo
            ncNew = ncOld + nc
            ! write(*,*) ncNew
            ! write(*,*) clusterList(ic,:)
        enddo
        ncTotal = ncNew
        ! write(*,*)
        ! write(*,*) 'i = ', i,'|', clusterList(ic,:)
    enddo
    k = 0
    nsub = 0.0
    do i = 1, icTotal
        nc = 0
        do j = 1, N
            if( clusterList(i,j) /= 0 )then
                nc = nc + 1
            endif
        enddo
        nsub = real(nc)**2/real(N) + nsub
        ! write(210+nRun,'(I3)', advance='no') nc
        k = k + nc
    enddo
    ! write(210+nRun,*) ''

    write(230+nRun,'(F7.2)', advance='no') nsub
    write(230+nRun,'(I7)', advance='no') tMCS-1
    write(230+nRun,*) ''
    if( k /= N )then
        write(*,*)'clusterSize: ERROR - you done goofed! k/=N'
    endif

end subroutine wrtClstrSize


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


! outputs x to fort.175
subroutine wrtXR( N, x, speciesR, tstep)
    ! rlxArea = preferred area of one cell
    ! L = number of lattice sites along one dimension
    ! N = total number of cells
    ! sigma = array of cell labels
    ! x = array of all the cells lattice sites
    implicit none
    integer,  intent(in) :: N, tstep
    integer,  intent(in), dimension(:,:,:) :: x
    real(b8), intent(in), dimension(:)     :: speciesR
    integer :: i, j

    ! write out x array
    do i = 1, N
        j = 1
        do while( x(i,j,1) /= 0 )
            write(175,*) x(i,j,1), x(i,j,2), speciesR(i), tstep-1
            j = j + 1
        enddo
    enddo

end subroutine wrtXR


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
