program main
    use func
    implicit none
    integer:: xgrid, ygrid
    real(8):: Tin, fuel
    integer:: dt
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

end program main
