program main
    use func
    implicit none
    integer:: xgrid, ygrid
    real(8):: Tin, fuel
    real(8):: dt, t
    integer:: max_gen, gen, print_gen
    real(8):: xlen, ylen
    real(8):: size, Thot
    real(8):: d
    character(16):: filename
    integer, parameter:: lout=125, lvtk=126
    real(8), allocatable, dimension(:,:):: U, V, Unew, Vnew
    real(8), allocatable, dimension(:)  :: X, Y

    xgrid = 500
    ygrid = xgrid
    Tin  = 1.39926
    fuel = 0.93781
    dt = 0.000001
    max_gen = 10000
    print_gen = 100
    gen = 0
    t = gen*dt
    xlen = 10.
    ylen = xlen
    size = 0.25
    Thot = Tin*3.
    d = 4.5
    filename = 'out.dat'

    allocate( X(xgrid) )
    allocate( Y(ygrid) )
    allocate( U(xgrid, ygrid) )
    allocate( V(xgrid, ygrid) )
    allocate( Unew(xgrid, ygrid) )
    allocate( Vnew(xgrid, ygrid) )

    X = linspace(-xlen/2, xlen/2, xgrid)
    Y = linspace(-ylen/2, ylen/2, ygrid)
    call initialize_uv(U, V, Tin, fuel)
    call make_hotspot(U, X, Y, size, Thot)
    !open(lout, file=filename, form='unformatted')
    !write(lout) t
    !write(lout) U
    !write(lout) V
    !write(lout) X
    !write(lout) Y
    !close(lout)
    call out_vtk(gen, U, X, Y, lvtk)

    do while ( gen <= max_gen )
        gen = gen + 1
        t = gen*dt

        call timestep(U, V, X, Y, d, dt, fuel, Unew, Vnew)
        call copy_2darr(Unew, U)
        call copy_2darr(Vnew, V)
        if ( mod(gen, print_gen) == 0 ) then
            write(*,*) "step =>",gen
            !call output(U, V, X, Y, t, filename, lout)
            call out_vtk(gen, U, X, Y, lvtk)
        end if

    end do

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

subroutine timestep(U, V, X, Y, d, dt, fuel, Unew, Vnew)
    real(8):: U(:,:), V(:,:), Unew(:,:), Vnew(:,:), X(:), Y(:)
    real(8):: dt
    real(8):: dx, dy
    real(8):: fuel
    real(8):: d
    integer:: i,j

    !.... Inner computation
    do i = 2, ubound(U,1)-1
        do j = 2, ubound(U,2)-1
            dx = ( X(i+1)-X(i-1) )/2.
            dy = ( Y(i+1)-Y(i-1) )/2.
            Unew(i, j) = U(i,j) + dt&
&                       *(laplacian(U(i-1,j), U(i,j), U(i+1,j), U(i,j-1), U(i,j+1), dx, dy)&
&                      + f(U(i,j), V(i,j)) )
            Vnew(i, j) = V(i,j) + dt&
&                      *(d*laplacian(V(i-1,j), V(i,j), V(i+1,j), V(i,j-1), V(i,j+1), dx, dy)&
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
end subroutine timestep

subroutine output(U, V, X, Y, t, filename, lout)
    real(8):: U(:,:), V(:,:), X(:), Y(:)
    real(8):: t
    character(16):: filename
    integer, intent(in):: lout

    open(lout, file=filename, form='unformatted', position='append')

    write(lout) t
    write(lout) U
    write(lout) V
    write(lout) X
    write(lout) Y
    close(lout)

    return
end subroutine output

subroutine out_vtk(step, U, X, Y, lvtk)
    integer, intent(in):: step
    integer, intent(in):: lvtk
    real(8), intent(in):: U(:,:)
    real(8), intent(in):: X(:), Y(:)
    character(40):: filename
    integer:: i,j

    write(filename, "(a,i5.5,a)") "data/temp",int(step),".vtk"
    open(lvtk, file=filename)
    write(lvtk, "('# vtk DataFile Version 3.0')")
    write(lvtk, "('test')")
    write(lvtk, "('ASCII ')")

    !.... Axis section
    write(lvtk, "('DATASET STRUCTURED_GRID')")
    write(lvtk, "('DIMENSIONS ', 3(1x,i5))") ubound(X,1), ubound(X,1), 1
    write(lvtk, "('POINTS ',i9,' float')") ubound(X,1)*ubound(Y,1)
    do i=1, ubound(X, 1)
        do j=1, ubound(Y, 1)
            write(lvtk, *) X(i), Y(j), 0.0
        end do
    end do

    !.... Data section
    write(lvtk, "('POINT_DATA ', i9)") ubound(U, 1) * ubound(U, 2)
    write(lvtk, "('SCALARS temp float')")
    write(lvtk, "('LOOKUP_TABLE default')")
    do i=1, ubound(U, 1)
        do j=1, ubound(U, 2)
            write(lvtk,*) U(i,j)
        end do
    end do
    close(lvtk)

    return
end subroutine out_vtk

end program main
