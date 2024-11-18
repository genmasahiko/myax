!This code is for vibration analysis of H2O molecule.

program disp
        implicit none

        integer, parameter :: dp = kind(0.0d0)
        character(100) :: input, output, new_input, new_input_dir
        character(2), allocatable :: asym(:)
        real(dp), allocatable :: apos(:,:)
        real(dp) :: dx 
        integer :: i, j, k, nat
        logical :: lfcp, tprnfor
        character(100) :: dum1, dum2

        write(*,*) 'relax input & output file name?'
        read(*,*) input, output

        nat = get_nat(input)
        
        allocate( asym(nat), apos(nat,3) )
        call read_output(output, nat, asym, apos)

        lfcp = .false. ; tprnfor = .false.
        call check_fcp_tprnfor(input, lfcp, tprnfor)

        write(*,*) lfcp, tprnfor

        do i = 1, nat
                write(*, '(a, 3f20.10)') asym(i), ( apos(i, j), j = 1, 3 )
        end do

        write(*,*) 'H2O shift for finite difference?'
        read(*,*) dx

        ! making 'vib' directory

        call system('mkdir vib')

        ! plus diff
        k = 0
        do i = nat - 2, nat
                do j = 1, 3
                        k = k + 1
                        write(new_input_dir, *) k
                        new_input_dir = 'vib/disp_'//trim(adjustl(new_input_dir))//'+'
                        new_input = trim(adjustl(new_input_dir))//'/in'

                        call system('mkdir ' //new_input_dir)

                        apos(i, j) = apos(i, j) + dx

                        call rewrite_input(input, new_input, asym, apos, nat, lfcp, tprnfor)

                        apos(i, j) = apos(i, j) - dx
                end do
        end do

        k = 0
        do i = nat - 2, nat
                do j = 1, 3
                        k = k + 1
                        write(new_input_dir, *) k
                        new_input_dir = 'vib/disp_'//trim(adjustl(new_input_dir))//'-'
                        new_input = trim(adjustl(new_input_dir))//'/in'

                        call system('mkdir ' //new_input_dir)

                        apos(i, j) = apos(i, j) - dx

                        call rewrite_input(input, new_input, asym, apos, nat, lfcp, tprnfor)

                        apos(i, j) = apos(i, j) + dx
                end do
        end do

contains

        function count_line(filename)
                implicit none

                character(100), intent(in) :: filename
                integer :: count_line, numline

                open(10, file=trim(filename), status="old")

                numline = 0

                do
                        read(10, *, end=100)
                        numline = numline + 1
                end do
100 close(10)

                close(10)

                count_line = numline
        end function

        subroutine read_output(filename, nat, atom_symbol, atom_position)
                implicit none

                integer, parameter :: dp = kind(0.0d0)
                integer :: i, j, nat
                character(100) :: filename
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

        subroutine rewrite_input(relax_input, new_input, atom_symbol, atom_position, nat, lfcp, tprnfor)
                implicit none

                character(100) :: relax_input, new_input, line
                character(2) :: atom_symbol(nat)
                real(dp) :: atom_position(nat,3)
                character(5) :: direction
                integer :: i, j, nat
                logical :: lfcp, tprnfor

                open(10, file=relax_input, status='old')
                open(11, file=new_input, status='replace')

                do
                        read(10, '(a)', end=100) line

                        if( index(line, 'control') > 0 .and. (.not. (tprnfor)) ) then
                                write(11, '(a)') trim(line)
                                write(11, '(a)') " tprnfor = .true."

                        else if ( index(line, 'ATOMIC_POSITIONS') > 0 ) then
                                if ( .not. lfcp ) then
                                        write(11, '(a)') trim(line)
                                        do i = 1, nat
                                                write(11, '(a, 3f20.10)') atom_symbol(i), &
                                                                ( atom_position(i, j), j = 1, 3 )
                                        end do
                                        exit
                                else
                                        write(11, '(a)') trim(line)
                                        do i = 1, nat
                                                write(11, '(a, 3f20.10, 3i4)') atom_symbol(i), &
                                                        ( atom_position(i, j), j = 1, 3 ), 0, 0, 0
                                        end do
                                        exit
                                end if

                        else if ( index(line, 'relax') > 0 .and. (.not. lfcp)) then
                                write(11, '(a)') " calculation = 'scf'"

                        else
                                write(11, '(a)') trim(line)
                        end if

                end do
100 close(10); close(11)

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

        subroutine check_fcp_tprnfor(filename, lfcp, tprnfor)
                implicit none

                character(100) :: filename, line
                integer :: ios
                logical :: lfcp, tprnfor

                open(10, file=filename, status='old')

                do
                        read(10, '(a)', iostat=ios) line
                        if ( ios /= 0 ) exit

                        if ( index(line, 'lfcp') > 0 ) then
                                lfcp = .true.
                        else if ( index(line, 'tprnfor') > 0 ) then
                                tprnfor = .true.
                        end if
                end do

                close(10)
        end subroutine

end program

