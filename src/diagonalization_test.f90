program main
        Implicit None

        integer, parameter :: dp=kind(0.0d0)
        real(dp), allocatable :: hes(:,:)
        integer :: i, j
        character(30) :: filename

        ! only for dgeev
        integer :: info
        real(dp), allocatable :: wr(:), wi(:), vl(:,:), vr(:,:), work(:)


        ! read a file containing Hessian
        write(*,*) "file name containing Hessian?"
        read(*,'(a)') filename

        open(unit=10, file=trim(filename), status='old', action='read')

        allocate( hes(9,9) )

        do i = 1, 9
                read(10,*) ( hes(i,j), j = 1, 9)
        end do

        close(10)

        write(*,*) "this matrix was read"
        do i = 1, 9
                write(*, '(9e20.8)') ( hes(i,j), j = 1, 9)
        end do

        allocate( wr(9), wi(9), vl(9,9), vr(9,9), work(8*9) )

        call dgeev('V', 'V', 9, hes, 9, wr, wi, vl, 9, vr, 9, work, 8*9, info)

        do i = 1, 9
                write(*,'(2f20.10)') wr(i), wi(i)
        end do

end program
