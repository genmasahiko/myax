program vibration
    implicit none
    integer, parameter :: dp = kind(0.0d0)

    real(dp), allocatable :: &
        & f_plus(:,:), f_minus(:,:), force(:,:,:,:),  &
        & hes(:,:), mass_hes(:,:), mass(:,:), &
        & omega(:), energy(:), wavenum(:), omega2(:), eig_vec(:,:), &
        & force_eq(:,:), force_average(:)

    integer :: i, j, k, nat, ios, l, pm
    character(100) :: filename, line, inp, input, output
    real(dp) :: dx, zpe
    real(dp) :: eigen_re(9), eigen_im(9)

    character(100), allocatable :: atom_symbol(:)
    real(dp), allocatable :: atom_position(:,:)

    real(dp) :: surf_position
    integer :: nat_surf


    ! Constants for unit conversion
    ! From UnitConverters.net
    real(dp), parameter :: Ang2Bohr = 1.8897259886d0 ! bohr/ang
    real(dp), parameter :: Amu2Me = 1822.8885300626d0 ! me/amu
    real(dp), parameter :: auTime2S = 2.4188843265857d-17 ! s/auTime, from Wiki
    real(dp), parameter :: pi = 3.14159265358979323846
    real(dp), parameter :: hbar = 6.582119514e-16 ! eV s
    real(dp), parameter :: c = 299792458 ! m/s

    inp = 'vib/disp_1+/in'
    nat = get_nat(inp)

    write(*,*) 'shift?'
    read(*,*) dx

    allocate( f_plus(9, 9), f_minus(9, 9) ) 
    allocate( force(nat, 3, 9, 2), force_average(3), force_eq(nat,3) )

    f_plus = 0
    f_minus = 0
    force = 0
    force_eq = 0

    ! for equilibrium position
    filename = 'out'
    call get_force( filename, nat, force_eq )

    do i = 1, 9
        do pm = 1, 2

            ! set filename 
            write(filename, *) i
            if ( pm == 1 ) then
                filename = 'vib/disp_'//trim(adjustl(filename))//'+'
            else if ( pm == 2 ) then
                filename = 'vib/disp_'//trim(adjustl(filename))//'-'
            end if
            input = trim(adjustl(filename))//'/out'

            call get_force( input, nat, force(:, :, i, pm) )
            !  call remove_drift_force( nat, force(:, :, i, pm), force_eq )

            if ( pm == 1 ) then
                l = 0
                do j = 1, 3
                    do k = 1, 3
                        l = l + 1
                        f_plus(i, l) = force( nat - 2 + j - 1, k, i, pm )
                    end do
                end do
            else if ( pm == 2 ) then
                l = 0
                do j = 1, 3
                    do k = 1, 3
                        l = l + 1
                        f_minus(i, l) = force( nat - 2 + j - 1, k, i, pm )
                    end do
                end do
            end if

        end do
    end do

    allocate( hes(9,9) )
    do i = 1, 9
        do j = 1, 9
            hes(i, j) = ( f_minus(i, j) - f_plus(i, j) ) / (2.0_dp*dx)
        end do
    end do

    deallocate( f_plus, f_minus )

    allocate( mass(9,9) )
    call get_mass_matrix(inp, nat, mass)

    allocate( mass_hes(9,9) )

    do i = 1, 9
        do j = 1, 9
            mass_hes(i, j) = hes(i, j) / mass(i, j)
        end do
    end do

    deallocate( mass, hes )

    allocate( eig_vec(9,9) )

    call diag(mass_hes, eigen_re, eigen_im, eig_vec)

    deallocate(mass_hes)

    write(*,*) '-----eigen value (Re, Im) of mass weighted hessian----- '
    do i = 1, 9
           write(*, '(2f20.10)') eigen_re(i), eigen_im(i)
    end do

    do i = 1, 9
        if ( dabs(eigen_re(i) ) < 1.0d-10 ) eigen_re(i) = 0.0_dp
    end do

    allocate( omega2(9), omega(9), energy(9), wavenum(9) )

    omega2 = 0.0_dp
    omega = 0.0_dp
    energy = 0.0_dp

    do i = 1, 9
        ! In Ry a.u., it must be multiplied by 2 because me/2 = 1.
        omega2(i) = eigen_re(i) / Ang2Bohr / ( Amu2Me * 2.0_dp )
        if (omega2(i) >= 0 ) then
            omega(i) = dsqrt(omega2(i)) / auTime2S
        else
            omega(i) = dsqrt(-omega2(i)) / auTime2S
        end if
        wavenum(i) = omega(i) * 1.0d-2 / (2*pi*c)
        energy(i) = hbar * omega(i)
    end do

    zpe = 0.0
    do i = 1, 9
        if ( omega2(i) >= 0 ) then
            zpe = zpe + energy(i) * 0.5_dp
        end if
    end do

    !deallocate( eigen_re, eigen_im )

    allocate( atom_symbol(nat), atom_position(nat, 3) )

    call read_output('out' , nat, atom_symbol, atom_position)

    ! calculate slab surface position for anime.x
    surf_position = 0.0_dp
    nat_surf = (nat - 3) / 4
    do i = nat - 3 - nat_surf + 1, nat - 3
        surf_position = surf_position + atom_position(i, 3)
    end do

    surf_position = surf_position / nat_surf

    output = 'vib/out.dat'

    open(10, file=output, status='replace')

    write(*,*) '----- frequency (Hz) & wave number (cm-1) & energy (eV) of H2O vibration -----'
    write(10,*) '----- frequency (Hz) & wave number (cm-1) & energy (eV) of H2O vibration -----'
    do i = 1, 9
        if (omega2(i) >= 0 ) then
            write(*, '(3f30.10)') omega(i), wavenum(i), energy(i)
            write(10, '(4f30.10)') omega(i), wavenum(i), energy(i)
        else
            write(*, '(f30.10,a, f29.10, a, f29.10, a)') omega(i), 'i', wavenum(i), 'i', energy(i), 'i'
            write(10, '(f30.10,a, f29.10, a, f29.10, a)') omega(i), 'i', wavenum(i), 'i', energy(i), 'i'
        end if
    end do

    write(*,*) '----------'
    write(10,*) '----------'

    write(*, '(a, f10.5, a6)') 'zero point energy = ', zpe, '(eV)'
    write(10, '(a, f10.5, a6)') 'zero point energy = ', zpe, '(eV)'

    write(*,*) '----- eigen vector -----'
    write(10,*) '----- eigen vector -----'
    do i = 1, 9
        write(*, '(9f16.10)') ( eig_vec(i, j), j = 1, 9 )
        write(10, '(9f16.10)') ( eig_vec(i, j), j = 1, 9 )
    end do

    write(*,*) '----- atom position -----'
    write(10,*) '----- atom position -----'
    do i = nat-2, nat
        write(*, '(3f16.10)') ( atom_position(i, j), j = 1, 3 )
        write(10, '(3f16.10)') ( atom_position(i, j), j = 1, 3 )
    end do

    write(*, '(a, f10.5)') 'surface position =', surf_position
    write(10, '(a, f10.5)') 'surface position =', surf_position

    close(10)

contains

    subroutine get_force(filename, nat, force)
        implicit none

        integer, parameter :: dp = kind(0.0d0)
        character(100) :: filename, line
        real(dp) :: force(nat, 3)
        integer :: ios, nat, i, j

        open(10, file=filename, status='old')

        do
            read(10, '(a)', iostat=ios) line
            if ( index(line, 'Forces acting on atoms') > 0 ) exit
        end do

        read(10, '(a)') ! read empty line

        do i = 1, nat
            read(10, '(a)') line
            line = line( index(line, '=') + 1:)
            line = trim(adjustl(line))
            read(line, *) ( force(i, j), j = 1, 3 )
        end do

        close(10)
    end subroutine

    subroutine remove_drift_force(nat, force, force_eq)
        implicit none

        integer, parameter :: dp = kind(0.0d0)
        real(dp), intent(inout) :: force(nat, 3)
        real(dp), intent(in) :: force_eq(nat, 3)
        real(dp) :: force_average(3)
        integer :: nat, i, j

        force_average = 0.0_dp

        ! substract equilibrium force
        ! NOTE: The necessity of this step is not clear.
        do i = 1, nat
            do j = 1, 3
                force(i, j) = force(i, j) - force_eq(i, j)
            end do
        end do

        ! calculate drift force as average of forces acting on H2O
        force_average = 0.0_dp
        do i = 1, 3
            do j = nat-2, nat
                force_average(i) = force_average(i) + force(j, i)
            end do
            force_average(i) = force_average(i) / 3.0_dp
        end do

        ! substract drift force
        do i = 1, 3
            do j = 1, nat
                force(j, i) = force(j, i) - force_average(i)
            end do
        end do

    end subroutine

    subroutine get_mass_matrix(filename, nat, mass)
        implicit none

        integer, parameter :: dp = kind(0.0d0)
        character(100) :: filename, line
        real(dp) :: mass(9,9)
        integer :: ios, nat, i, j, k
        character(100) :: symbol(nat)
        real(dp) :: m(9)


        ! atomic weight from CIAAW of IUPAC
        real(dp), parameter :: mass_H = 1.008, mass_O = 15.999


        open(10, file=filename, status='old')

        do
            read(10, '(a)', iostat=ios) line
            if ( ios == -1 ) exit

            if ( index(line, 'ATOMIC_POSITIONS') > 0 ) then
                do i = 1, nat
                    read(10, *) symbol(i)
                end do
                exit
            end if
        end do

        m = 0

        j = 0
        do i = nat - 2, nat
            if ( trim(symbol(i)) == 'H' ) then
                do k = 1, 3
                    j = j + 1
                    m(j) = dsqrt(mass_H) 
                end do
            else if ( trim(symbol(i)) == 'O' ) then
                do k = 1, 3
                    j = j + 1
                    m(j) = dsqrt(mass_O)
                end do
            end if
        end do

        close(10)

        do i = 1, 9
            do j = 1, 9
                mass(i, j) = m(i) * m(j)
            end do
        end do
    end subroutine

    function get_nat(filename)
        implicit none

        character(100) :: filename, line
        character(5) :: dum1, dum2
        integer :: get_nat

        open(10, file=filename, status='old')

        do
            read(10, '(a)') line
            if ( index(line, 'nat') > 0 ) then
                read(line, *) dum1, dum2, get_nat
                exit
            end if
        end do

        close(10)
    end function

    subroutine diag(matrix, eigen_re, eigen_im, vr)
        implicit none

        integer, parameter :: dp=kind(0.0d0)
        real(dp) :: matrix(9, 9)

        ! only for dgeev
        integer :: info
        real(dp), allocatable :: vl(:,:), work(:)
        real(dp) :: vr(9,9)
        real(dp) :: eigen_re(9), eigen_im(9)

        allocate( vl(9,9), work(8*9) )
        call dgeev('V', 'V', 9, matrix, 9, eigen_re, eigen_im, vl, 9, vr, 9, work, 8*9, info)
        deallocate( vl, work )
    end subroutine

    subroutine read_output(filename, nat, atom_symbol, atom_position)
        implicit none

        integer, parameter :: dp = kind(0.0d0)
        integer :: i, j, nat
        character(3) :: filename
        character(100) :: line
        character(2) :: atom_symbol(nat)
        real(dp) :: atom_position(nat,3)

        open(10, file=trim(filename), status='old')

        do
            read(10, '(a)', end=100) line
            if ( index(line, 'ATOMIC_POSITIONS') > 0 ) then
                do i = 1, nat
                    read(10, *) atom_symbol(i), ( atom_position(i, j), j = 1, 3 )
                end do
            end if
        end do
        100 close(10)
    end subroutine

end program
