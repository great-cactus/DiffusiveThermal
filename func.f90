module func
implicit none
public
    real(8), save, private:: eps, a, d, h
    logical, save, public :: isConstRead = .false.

contains

logical function const_read()
    const_read = isConstRead
end function

!##
!.... Here are the functions used for computation
real(8) function k(u, j, j)
    real(8):: u
    integer:: i, j

    k = exp( u(i,j) / (1 + u(j,j)/c) )
end function

real(8) function f(u, v, i, j)
    real(8):: u, v
    integer:: i, j

    f = 1./eps * ( -a*u(i,j) + k( u(i,j) )*v(i,j) )
end function

real(8) function g(u, v, i, j)
    real(8):: u, v
    integer:: i, j

    g = h*( 1-v(i,j) ) -k( u(i,j) ) * v(i,j)
end function

real(8) laplacian(f, i, j, di, dj)
    real(8):: f
    real(8):: di, dj
    integer:: i, j

    laplacian = ( f(i+1,j)-2*f(i,j)+f(i-1,j) )/(di*di)&
    &         + ( f(i,j+1)-2*f(i,j)+f(j,i-1) )/(dj*dj)
end function
!.... end
!##




end module func
