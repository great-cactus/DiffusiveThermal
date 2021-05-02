program main
    use func
    implicit none
    integer:: ngrid=11
    integer:: dt
    real(8):: xlen=10., ylen=10.
    real(8), allocatable, dimension(:,:):: U, V
    real(8), allocatable, dimension(:)  :: X, Y

    allocate( X(ngrid) )
    allocate( Y(ngrid) )
    allocate( U(ngrid, ngrid) )
    allocate( V(ngrid, ngrid) )

    X = linspace(-xlen/2, xlen/2, ngrid)
    Y = linspace(-ylen/2, ylen/2, ngrid)
    print *,X(:)

end program main
