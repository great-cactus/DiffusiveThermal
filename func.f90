module func
    implicit none
    public
    real(8), save, private:: eps = 0.0125, a = 2., d = 4.5, h = 45., c = 5.
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
    integer:: n, i, ierr
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
