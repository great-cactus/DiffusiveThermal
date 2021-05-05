module func
    implicit none
    public
    real(8), save, private:: eps = 0.0005, a = 2., h = 45., c = 5.
    logical, save, public :: isConstRead = .false.

contains

logical function const_read()
    const_read = isConstRead
    return
end function const_read

!##
!.... Here are the functions used for computation
real(8) function k(u_ij)
    real(8):: u_ij

    k = exp( u_ij / (1 + u_ij/c) )
    return
end function k

real(8) function f(u_ij, v_ij)
    real(8):: u_ij, v_ij

    f = 1./eps * ( -a*u_ij + k( u_ij )*v_ij )
    return
end function f

real(8) function g(u_ij, v_ij)
    real(8):: u_ij, v_ij

    g = h*( 1-v_ij ) -k( u_ij ) * v_ij
    return
end function g

real(8) function laplacian(fi_bw, fij, fi_fw, fj_bw, fj_fw, di, dj)
    real(8):: fi_bw, fij, fi_fw, fj_bw, fj_fw
    real(8):: di, dj

    laplacian = ( fi_fw- 2*fij + fi_bw )/(di*di)&
    &         + ( fj_fw- 2*fij + fj_bw )/(dj*dj)
    return
end function laplacian
!.... end
!##

!##
!.... Here is the functions can be used for simple coding
function linspace(x0, x1, nd) result(x)
!.... Makes array of nx number of points
    real(8), intent(in):: x0, x1
    integer, optional, intent(in):: nd
    real(8), dimension(:), allocatable:: x
    integer:: n, i
    real(8):: dx

    if ( present(nd) ) then
        n = nd
        if (n < 1) then
            call stop_error("meshexp: nx > 1 required.")
        endif
    else
        n=100!    Default n
    endif

    allocate(x(n))
    if (n == 1) then
        x(1) = x1
        return
    endif

    dx = (x1-x0)/dble(n-1)
    do i=1, n
        x(i) = x0 + dx*dble(i-1)
    enddo

    return
end function linspace

subroutine copy_2darr(arr, newarr)
    real(8):: arr(:,:), newarr(:,:)
    integer:: i, j

    if (( ubound(arr, 1) /= ubound(newarr, 1) ) .or. ( ubound(arr, 2) /= ubound(newarr, 2) )) then
        call stop_error('arr and newarr has different size.')
    end if

    do i = 1, ubound(arr,1)
        do j = 1, ubound(arr,2)
            newarr(i, j) = arr(i, j)
        enddo
    enddo

    return
end subroutine copy_2darr

real(8) function rand_gen()
    integer, allocatable, dimension(:):: seed
    integer:: i, seedsize
    character( 8):: date
    character(10):: time
    character( 5):: zone
    integer      :: values(8)

    call random_seed(size=seedsize)
    allocate( seed(seedsize) )
    call random_seed( get=seed )
    call date_and_time( date, time, zone, values )
    seed(seedsize) = values(8)
    call random_seed(put=seed)

    call random_number(rand_gen)

    return
end function rand_gen

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
! stop_error() uses the statement "stop 1" which returns an exit code
! 1 and prints the messeage.

    character(len=*):: msg
    print *, msg
    stop 1
end subroutine stop_error
!.... end
!##

end module func
