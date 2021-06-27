program main
    use func
    implicit none
    real(8):: eps, a, h, c
    integer:: xgrid, ygrid
    real(8):: Tin, fuel
    real(8):: dt, t
    integer:: max_gen, gen, print_gen
    real(8):: xlen, ylen
    real(8):: size, Thot
    integer:: n_hot
    real(8):: d
    integer, parameter:: lout=125
    character(20)     :: out_file='DT.out'
    integer:: lrnd
    real(8), allocatable, dimension(:,:):: U, V, Unew, Vnew
    real(8), allocatable, dimension(:)  :: X, Y

    xgrid = 500
    ygrid = xgrid
    Tin  = 1.39926
    fuel = 0.93781
    dt = 0.000001
    max_gen =   100000
    print_gen = 1000
    gen = 0
    t = gen*dt
    xlen = 10.
    ylen = xlen
    n_hot = 3
    size = 0.2
    Thot = 10.
    lrnd = 1
    d = 2.
    a = 2
    h = 45
    c = 5
    eps = 0.001

    call read_inp(a, c, h, eps, d,&
                 &xlen, ylen, xgrid, ygrid, Tin, fuel,&
                 &max_gen, print_gen, dt,&
                 &n_hot, lrnd, size, Thot)

    allocate( X(xgrid) )
    allocate( Y(ygrid) )
    allocate( U(xgrid, ygrid) )
    allocate( V(xgrid, ygrid) )
    allocate( Unew(xgrid, ygrid) )
    allocate( Vnew(xgrid, ygrid) )

    X = linspace(-xlen/2, xlen/2, xgrid)
    Y = linspace(-ylen/2, ylen/2, ygrid)
    call initialize_uv(U, V, Tin, fuel)
    if (lrnd == 1) then
        call make_n_spots(U, X, Y, size, Thot, n_hot)
    else if (lrnd == 0) then
        call make_hotspot(U, X, Y, size, Thot)
    else
        call stop_error('Invarid lrand.')
    end if

    call out_vtk(gen, U, X, Y)
    do while ( gen <= max_gen )
        gen = gen + 1
        t = gen*dt

        call timestep(U, V, X, Y, d, dt, Unew, Vnew)
        call copy_2darr(Unew, U)
        call copy_2darr(Vnew, V)
        if ( mod(gen, print_gen) == 0 ) then
            open (lout, file=out_file, form='formatted',&
                 &position='append', action='write', status='old')

                write(lout,*) "step =>",gen

            close(lout)
            call out_vtk(gen, U, X, Y)
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
    i_max = int( ubound(U, 1)*0.8 )
    i_min = int( ubound(U, 1)*0.2 )
    j_max = int( ubound(U, 2)*0.8 )
    j_min = int( ubound(U, 2)*0.2 )

    do k=1, n
        centor_idx(k, 1) = int( rand_gen() * (i_max - i_min) + i_min )
        centor_idx(k, 2) = int( rand_gen() * (j_max - j_min) + j_min )
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

subroutine timestep(U, V, X, Y, d, dt, Unew, Vnew)
    real(8):: U(:,:), V(:,:), Unew(:,:), Vnew(:,:), X(:), Y(:)
    real(8):: dt
    real(8):: dx, dy
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

subroutine out_vtk(step, U, X, Y)
    implicit none
    integer, intent(in):: step
    real(8), intent(in):: U(:,:)
    real(8), intent(in):: X(:), Y(:)
    integer, parameter:: lvtk=126
    character(40):: filename
    integer:: i,j

    write(filename, "(a,i8.8,a)") "data/temp",int(step),".vtk"
    open(lvtk, file=filename, action='write',&
        &form='formatted', position='rewind')

        !.... Initialize vtk
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

subroutine read_inp(a, c, h, eps, d,&
                   &xlen, ylen, xgrid, ygrid, Tin, fuel,&
                   &max_gen, print_gen, dt,&
                   &n_hot, lrnd, size, Thot)
    implicit none
    integer:: line, ierr, pos
    character(100):: buffer, keyword, err_msg
    real(8):: eps, a, h, c
    integer:: xgrid, ygrid
    real(8):: Tin, fuel
    real(8):: dt
    integer:: max_gen, print_gen
    real(8):: xlen, ylen
    real(8):: size, Thot
    integer:: n_hot
    real(8):: d
    integer:: lrnd
    integer, parameter:: linp=124,          lout=125
    character(20)     :: inp_file='DT.inp', out_file='DT.out'

    open (lout, file=out_file, form='formatted',&
         &position="rewind", action='write')

        write(lout,"(A)")"*****************************************************"
        write(lout,"(A)")""
        write(lout,"(A)")"Starting 2D diffusive thermal computations."
        write(lout,"(A)")"This code was made by akira tsunoda."
        write(lout,"(A)")""
        write(lout,"(A)")"*****************************************************"

    close (lout)

    open(linp, file=inp_file, form='formatted',&
        &position='rewind', action='read', status='old')
    open(lout, file=out_file, form='formatted',&
        &position="append", action='write', status='old')

        write(lout, "(A)") ""
        write(lout, "(A)") "Start reading input file.........."
        write(lout, "(A)") ""

        ierr = 0
        do while (ierr == 0)
            read(linp, '(A)', iostat=ierr) buffer
            write(lout, "(A)") buffer
            if (ierr == 0) then
                line = line + 1
                pos = scan(buffer, '    ')
                keyword = buffer(1:pos)
                buffer = buffer(pos+1:)
                select case (keyword)
                    case ('#')!    Comment
                        continue
                    !##
                    !##    Constants for functions
                    !##
                    !   Constant a
                    case ('CONA')
                        read(buffer, *, iostat=ierr) a
                    !   Constant c
                    case ('CONC')
                        read(buffer, *, iostat=ierr) c
                    !   Constant h
                    case ('CONH')
                        read(buffer, *, iostat=ierr) h
                    !   Constant epsilon
                    case ('CEPS')
                        read(buffer, *, iostat=ierr) eps
                    !   Constant d
                    case ('COND')
                        read(buffer, *, iostat=ierr) d
                    !##
                    !##    Domain parameters
                    !##
                    !   X computation domain size
                    case ('XLEN')
                        read(buffer, *, iostat=ierr) xlen
                    !   Y computation domain size
                    case ('YLEN')
                        read(buffer, *, iostat=ierr) ylen
                    !   X point numbers
                    case ('XNUM')
                        read(buffer, *, iostat=ierr) xgrid
                    !   Y point numbers
                    case ('YNUM')
                        read(buffer, *, iostat=ierr) ygrid
                    !   Inlet U
                    case ('UINT')
                        read(buffer, *, iostat=ierr) Tin
                    !   Inlet V
                    case ('VINT')
                        read(buffer, *, iostat=ierr) fuel
                    !##
                    !##    Timestep parameters
                    !##
                    !   Max generation
                    case ('GENM')
                        read(buffer, *, iostat=ierr) max_gen
                    !   Print generation
                    case ('GENP')
                        read(buffer, *, iostat=ierr) print_gen
                    !   computation timestep
                    case ('DT')
                        read(buffer, *, iostat=ierr) dt
                    !##
                    !##    Initial conditions
                    !##
                    !   Initial hotspots number
                    case ('NHOT')
                        read(buffer, *, iostat=ierr) n_hot
                    !   0: hotspot at centor, 1: random hotspot
                    case ('LRND')
                        read(buffer, *, iostat=ierr) lrnd
                    !   size of hotspot
                    case ('SHOT')
                        read(buffer, *, iostat=ierr) size
                    !   magnitude of hotspot
                    case ('THOT')
                        read(buffer, *, iostat=ierr) Thot

                    case default
                        write(err_msg, "( a,i2 )") 'invalid keyword detected at line', line
                        call stop_error(err_msg)
                end select
            end if
        end do

        write(lout, "(A)") ""
        write(lout, "(A)") "Finish reading input file successfully."
        write(lout, "(A)") ""
    close(lout)
    close(linp)

    call read_consts(eps, a, h, c)

    return
end subroutine read_inp

end program main
