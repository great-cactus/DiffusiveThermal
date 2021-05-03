program main
    use func
    implicit none
    integer:: xgrid, ygrid
    real(8):: Tin, fuel
    real(8):: dt, t, t_end
    real(8):: xlen, ylen
    real(8):: size, Thot
    real(8), allocatable, dimension(:,:):: U, V
    real(8), allocatable, dimension(:)  :: X, Y

    xgrid = 11
    ygrid = xgrid
    Tin  = 1.39926
    fuel = 0.93781
    dt = 0.00001
    xlen = 10.
    ylen = xlen
    size = 2
    Thot = 2

    allocate( X(xgrid) )
    allocate( Y(ygrid) )
    allocate( U(xgrid, ygrid) )
    allocate( V(xgrid, ygrid) )

    X = linspace(-xlen/2, xlen/2, xgrid)
    Y = linspace(-ylen/2, ylen/2, ygrid)
    call initialize_uv(U, V, Tin, fuel)
    call make_hotspot(U, X, Y, size, Thot)

    stop
contains

subroutine initialize_uv(U, V, Tin, fuel)
    real(8):: U(:,:), V(:,:)
    real(8):: Tin, fuel
    integer:: i,j

    do i = 1, ubound(V,1)
        do j = 1, ubound(V,2)
            U(i, j) = Tin
            V(i, j) = fuel
        enddo
    enddo

    return
end subroutine initialize_uv

subroutine make_hotspot(U, X, Y, size, Thot)
    real(8):: U(:,:), X(:), Y(:)
    real(8):: size, Thot
    integer::i, j

    do i = 1, ubound(U,1)
        do j = 1, ubound(U,2)
            if( ( abs(X(i)) <= size/2. ) .and. ( abs(Y(j)) <= size/2. ) ) then
                U(i, j) = Thot
            endif
        enddo
    enddo

    return
end subroutine make_hotspot

real(8) function timestep(U, V, X, Y, dt, fuel) result(Unew, Vnew)
    real(8):: U(:,:), V(:,:), X(:), Y(:)
    real(8):: dt
    real(8):: dx, dy
    real(8):: fuel
    real(8), dimension(ubound(U,1),ubound(U,2)):: Unew
    real(8), dimension(ubound(V,1),ubound(V,2)):: Vnew

    !.... Inner computation
    do i = 2, ubound(U,1)-1
        do j = 2, ubound(U,2)-1
            dx = ( X(i+1)-X(i-1) )/2.
            dy = ( Y(i+1)-Y(i-1) )/2.
            Unew(i, j) = U(i,j) + dt&
&                        (laplacian(U(i-1,j), U(i,j), U(i+1,j), U(i,j-1), U(i,j+1), dx, dy)&
&                      + f(U(i,j), V(i,j)) )
            Vnew(i, j) = V(i,j) + dt&
&                        (laplacian(V(i-1,j), V(i,j), V(i+1,j), V(i,j-1), V(i,j+1), dx, dy)&
&                      + g(U(i,j), V(i,j)) )
        enddo
    enddo

    !.... Boundary computation
    do i = 1, ubound(U,1)
        Unew(i, 1)             = U(i, 2)
        Unew( i, ubound(U,2) ) = U( i, ubound(U,2)-1 )
        Vnew(i, 1)             = fuel
        Vnew( i, ubound(V,2) ) = fuel
    enddo
    do j = 1, ubound(U,2)
        Unew(1, j)             = U(1, j)
        Unew( ubound(U,1), j ) = U( ubound(U,1)-1, j )
        Vnew(1, j)             = fuel
        Vnew( ubound(U,1), j ) = fuel
    enddo

    return
end function timestep

end program main
