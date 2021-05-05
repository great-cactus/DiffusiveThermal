program main
    use func
    implicit none
    integer:: xgrid, ygrid
    real(8):: Tin, fuel
    real(8):: dt, t
    integer:: max_gen, gen, print_gen
    real(8):: xlen, ylen
    real(8):: size, Thot
    integer:: n_hot
    real(8):: d
    integer, parameter:: lvtk=126
    real(8), allocatable, dimension(:,:):: U, V, Unew, Vnew
    real(8), allocatable, dimension(:)  :: X, Y

    xgrid = 500
    ygrid = xgrid
    Tin  = 1.39926
    fuel = 0.93781
    dt = 0.000001
    max_gen =   1000000
    print_gen = 10000
    gen = 0
    t = gen*dt
    xlen = 10.
    ylen = xlen
    n_hot = 3
    size = 0.2
    Thot = Tin*3.
    d = 9.

    allocate( X(xgrid) )
    allocate( Y(ygrid) )
    allocate( U(xgrid, ygrid) )
    allocate( V(xgrid, ygrid) )
    allocate( Unew(xgrid, ygrid) )
    allocate( Vnew(xgrid, ygrid) )

    X = linspace(-xlen/2, xlen/2, xgrid)
    Y = linspace(-ylen/2, ylen/2, ygrid)
    call initialize_uv(U, V, Tin, fuel)
    !call make_hotspot(U, X, Y, size, Thot)
    call make_n_spots(U, X, Y, size, Thot, n_hot)
    call out_vtk(gen, U, X, Y, lvtk)

    do while ( gen <= max_gen )
        gen = gen + 1
        t = gen*dt

        call timestep(U, V, X, Y, d, dt, fuel, Unew, Vnew)
        call copy_2darr(Unew, U)
        call copy_2darr(Vnew, V)
        if ( mod(gen, print_gen) == 0 ) then
            write(*,*) "step =>",gen
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

subroutine make_n_spots(U, X, Y, size, Thot, n)
    real(8):: U(:,:), X(:), Y(:)
    real(8):: size, Thot
    integer:: n
    integer, dimension(n,2):: centor_idx
    integer:: i_max, i_min, j_max, j_min
    real(8):: x_centor, y_centor
    integer:: i, j, k

    !.... Define idxs to be centors
    i_max = ubound(U, 1)*0.8
    i_min = ubound(U, 1)*0.2
    j_max = ubound(U, 2)*0.8
    j_min = ubound(U, 2)*0.2

    do k=1, n
        centor_idx(k, 1) = rand_gen() * (i_max - i_min) + i_min
        centor_idx(k, 2) = rand_gen() * (j_max - j_min) + j_min
    end do

    !.... Make hot spot
    do k = 1, n
        x_centor = X( centor_idx(k, 1) )
        y_centor = Y( centor_idx(k, 2) )
        do i = 1, ubound(U,1)
            do j = 1, ubound(U,2)
                if( ( abs(X(i)-x_centor) <= size/2. )&
              &.and.( abs(Y(j)-y_centor) <= size/2. ) ) then
                    U(i, j) = Thot
                endif
            enddo
        enddo
    end do

    return
end subroutine make_n_spots

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
        Vnew(i, 1)             = V(i, 2)
        Vnew( i, ubound(V,2) ) = V( i, ubound(V,2)-1 )
    enddo
    do j = 1, ubound(U,2)
        Unew(1, j)             = U(2, j)
        Unew( ubound(U,1), j ) = U( ubound(U,1)-1, j )
        Vnew(1, j)             = V(2, j)
        Vnew( ubound(U,1), j ) = V( ubound(V,1)-1, j )
    enddo

    return
end subroutine timestep

subroutine out_vtk(step, U, X, Y, lvtk)
    integer, intent(in):: step
    integer, intent(in):: lvtk
    real(8), intent(in):: U(:,:)
    real(8), intent(in):: X(:), Y(:)
    character(40):: filename
    integer:: i,j

    write(filename, "(a,i8.8,a)") "data/temp",int(step),".vtk"
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
    write(lvtk, "('SCALARS U float')")
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
