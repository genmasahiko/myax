program anime
        implicit none

        integer, parameter :: dp = kind(0.0d0)
        real(dp), parameter :: pi = 3.14159265358979

        real(dp) :: dt, tmax
        integer :: tstep

        character(100) :: filename, vibout
        integer :: i, j, im(9)
        real(dp) :: y
        real(dp) :: freq(9), vec(9,9), r0(9), r1(9), r(9), surf

        vibout = 'out.dat'
        call read_vibout(vibout, freq, vec, r0, surf, im)

        tstep = 100

        do i = 1, 9
                if ( im(i) == 1 ) then
                        tmax = 2.0 * pi / freq(i)
                else
                        tmax = 2.0 * pi * 6.0 / freq(i)
                end if

                dt = tmax / tstep

                write(*,'(a, i2)') 'mode', i

                call time_evo(freq, vec, r0, r1, surf, r, tstep, dt, i, im(i) )
                call make_plt_xy(i, dt, tstep, r1, im(i))
                call make_plt_zx(i, dt, tstep, r1, im(i))

        end do

contains
        subroutine read_vibout(filename, freq, vec, r0, surf, im)
                implicit none

                integer, parameter :: dp = kind(0.0d0)

                character(100) :: filename, line
                integer :: i, im(9)
                real(dp) :: vec(9,9), r0(9), freq(9), surf

                open(10, file=filename, status='old')

                do
                        read(10, '(a)') line

                        if ( index(line, 'frequency') > 0 ) then
                                do i = 1, 9
                                        read(10, '(a)') line
                                        if ( index(line, 'i') > 0 ) then
                                                im(i) = 1
                                                line = line( 1: index(line, 'i')-1 ) 
                                                read(line, *) freq(i)
                                        else
                                                im(i) = 0
                                                read(line, *) freq(i)
                                        end if
                                end do

                        else if ( index(line, 'eigen vector') > 0 ) then
                                read(10, *) vec
                        
                        else if ( index(line, 'atom position') > 0 ) then
                                read(10, *) r0

                        else if ( index(line, 'surface position') > 0 ) then
                                if ( index(line, 'NaN') > 0 ) then
                                        surf = 0
                                        exit
                                else
                                        line = line( index(line, '=')+1 : )
                                        read(line, *) surf
                                        exit
                                end if
                        end if
                end do

                close(10)
        end subroutine

        subroutine time_evo(freq, vec, r0, r1, surf, r, tstep, dt, mode_num, ilabel)
                implicit none

                integer, parameter :: dp = kind(0.0d0)

                integer :: i, tstep, mode_num, ilabel
                real(dp) :: freq(9), vec(9,9), r0(9), r(9), surf, r1(9)
                real(dp) :: dt, time

                if ( access('./anime/', 'r') > 0 ) then
                        call system('mkdir anime')
                end if

                write(filename, *) mode_num
                filename = 'anime/'//trim(adjustl(filename))//'.dat' 


                ! r1 is the position shifted so that the x(O) = y(O) = 0, 
                ! and z(O) is distance between surface and O

                do i = 1, 9
                        if ( mod(i, 3) > 0 ) then
                                r1(i) = r0(i) - r0( mod(i - 1, 3) + 1 )
                        else
                                r1(i) = r0(i) - surf
                        end if
                end do


                open(10, file=filename, status='replace')

                time = 0.0

                if ( ilabel == 1 ) then
                        do i = 1, tstep
                                time = i * dt

                                do j = 1, 9
                                        if ( j <= 3 ) then
                                                r(j) = r1(j) + vec(mode_num, j) / dsqrt(15.999_dp) * dexp( - freq(mode_num) * time )
                                        else
                                                r(j) = r1(j) + vec(mode_num, j) / dsqrt(1.008_dp) * dexp( - freq(mode_num) * time )
                                        end if
                                end do

                                write(10, '(e20.10, 9f20.10)') time, ( r(j), j = 1, 9 )
                        end do
                else 
                        do i = 1, tstep
                                time = i * dt

                                do j = 1, 9
                                        if ( j <= 3 ) then
                                                r(j) = r1(j) + vec(mode_num, j) / dsqrt(15.999_dp) * dcos( freq(mode_num) * time )
                                        else
                                                r(j) = r1(j) + vec(mode_num, j) / dsqrt(1.008_dp) * dcos( freq(mode_num) * time )
                                        end if
                                end do

                                write(10, '(e20.10, 9f20.10)') time, ( r(j), j = 1, 9 )
                        end do
                end if

                close(10)
                        
        end subroutine

        subroutine make_plt_xy(mode_num, dt, tstep, r1, ilabel)
                implicit none

                integer, parameter :: dp = kind(0.0d0)

                character(100) :: pltfile, output, input
                integer :: tstep, mode_num, ilabel
                real(dp) :: dt
                real(dp) :: r1(9)
                character(100) :: zmin, zmax

                pltfile = 'anime/plot.plt'

                write(zmin, '(f5.1)') r1(3) - 2.0
                write(zmax, '(f5.1)') r1(3) + 2.0

                write(input, *) mode_num
                input = 'anime/'//trim(adjustl(input))//'.dat'

                write(output, *) mode_num
                output = 'anime/'//trim(adjustl(output))//'_xy.gif'
                
                open(10, file=pltfile, status='replace')

                write(10, *) 'set terminal gif animate delay 10 size 400,400'
                write(10, *) 'set output ', '"', trim(adjustl(output)), '"'
                write(10, *) 'set xlabel "x(Å)"; set ylabel "y(Å)" offset 2,0; set zlabel ""'
                write(10, *) 'set ztics'
                write(10, *) 'set format x "%.1f"'
                write(10, *) 'set format y "%.1f"'
                write(10, *) 'set xtics 1.0'
                write(10, *) 'set ytics 1.0'
                write(10, *) 'set xrange [-2:2]; set yrange [-2:2]; set zrange [', trim(adjustl(zmin)), ':', &
                        trim(adjustl(zmax)), ']'
                write(10, *) 'set view 0, 0'
                write(10, *) 'tstep = ', tstep
                write(10, *) 'dt = ', dt
                write(10, *) 'set parametric'
                write(10, *) 'set samples 18'
                write(10, *) 'set isosamples 18'
                write(10, *) 'set hidden3d'
                write(10, *) 'set view equal xyz'
                write(10, *) 'Fx(u,v)=sin(u)*cos(v)'
                write(10, *) 'Fy(u,v)=sin(u)*sin(v)'
                write(10, *) 'Fz(u,v)=cos(u)'
                write(10, *) "com1 = sprintf(", '"awk ', "'{print $2, $3, $4}' ", trim(adjustl(input)), '")'
                write(10, *) "com2 = sprintf(", '"awk ', "'{print $5, $6, $7}' ", trim(adjustl(input)), '")'
                write(10, *) "com3 = sprintf(", '"awk ', "'{print $8, $9, $10}' ", trim(adjustl(input)), '")'
                write(10, *) 'o = system(com1)'
                write(10, *) 'h1 = system(com2)'
                write(10, *) 'h2 = system(com3)'
                write(10, *) 'or = 0.74*0.4'
                write(10, *) 'hr = 0.46*0.4'
                write(10, *) 'do for [i=1:tstep] {'
                write(10, *) 'ox = real(word(o, 3*i-2))'
                write(10, *) 'oy = real(word(o, 3*i-1))'
                write(10, *) 'oz = real(word(o, 3*i))'
                write(10, *) 'h1x = real(word(h1, 3*i-2))'
                write(10, *) 'h1y = real(word(h1, 3*i-1))'
                write(10, *) 'h1z = real(word(h1, 3*i))'
                write(10, *) 'h2x = real(word(h2, 3*i-2))'
                write(10, *) 'h2y = real(word(h2, 3*i-1))'
                write(10, *) 'h2z = real(word(h2, 3*i))'
                if ( ilabel == 1 ) then
                        write(10, *) 'set title sprintf("Time = %.2e\n Imaginary frequency", i * dt)'
                else
                        write(10, *) 'set title sprintf("Time = %.2e", i * dt)'
                end if
                write(10, *) 'set urange [0:pi]'
                write(10, *) 'set vrange [0:2*pi]'
                write(10, *) 'splot or*Fx(u,v)+ox, or*Fy(u,v)+oy, or*Fz(u,v)+oz linetype rgbcolor "0xFE0300" notitle, \'
                write(10, *) 'hr*Fx(u,v)+h1x, hr*Fy(u,v)+h1y, hr*Fz(u,v)+h1z linetype rgbcolor "0xFFCCCC" notitle, \'
                write(10, *) 'hr*Fx(u,v)+h2x, hr*Fy(u,v)+h2y, hr*Fz(u,v)+h2z linetype rgbcolor "0xFF7577" notitle'
                write(10, *) '}'

                close(10)

                call system('gnuplot "anime/plot.plt"')

        end subroutine

        subroutine make_plt_zx(mode_num, dt, tstep, r1, ilabel)
                implicit none

                integer, parameter :: dp = kind(0.0d0)

                character(100) :: pltfile, output, input
                integer :: tstep, mode_num, ilabel
                real(dp) :: dt
                real(dp) :: r1(9)
                character(100) :: zmin, zmax

                pltfile = 'anime/plot.plt'

                write(zmin, '(f5.1)') r1(3) - 2.0
                write(zmax, '(f5.1)') r1(3) + 2.0

                write(input, *) mode_num
                input = 'anime/'//trim(adjustl(input))//'.dat'

                write(output, *) mode_num
                output = 'anime/'//trim(adjustl(output))//'_zx.gif'
                
                open(10, file=pltfile, status='replace')

                write(10, *) 'set terminal gif animate delay 10 size 400,400'
                write(10, *) 'set output ', '"', trim(adjustl(output)), '"'
                write(10, *) 'set label "x(Å)" at 0, 0, ', trim(adjustl(zmin)) ,'-3; set ylabel ""; set zlabel "z(Å)"'
                write(10, *) 'set ytics ("""")'
                write(10, *) 'set xrange [-2:2]; set yrange [-2:2]; set zrange [', trim(adjustl(zmin)), ':', &
                        trim(adjustl(zmax)), ']'
                write(10, *) 'set format x "%.1f"'
                write(10, *) 'set format z "%.1f"'
                write(10, *) 'set xtics 1.0'
                write(10, *) 'set ztics 1.0'
                write(10, *) 'set view 90, 0'
                write(10, *) 'tstep = ', tstep
                write(10, *) 'dt = ', dt
                write(10, *) 'set parametric'
                write(10, *) 'set samples 18'
                write(10, *) 'set isosamples 18'
                write(10, *) 'set hidden3d'
                write(10, *) 'set view equal xyz'
                write(10, *) 'Fx(u,v)=sin(u)*cos(v)'
                write(10, *) 'Fy(u,v)=sin(u)*sin(v)'
                write(10, *) 'Fz(u,v)=cos(u)'
                write(10, *) "com1 = sprintf(", '"awk ', "'{print $2, $3, $4}' ", trim(adjustl(input)), '")'
                write(10, *) "com2 = sprintf(", '"awk ', "'{print $5, $6, $7}' ", trim(adjustl(input)), '")'
                write(10, *) "com3 = sprintf(", '"awk ', "'{print $8, $9, $10}' ", trim(adjustl(input)), '")'
                write(10, *) 'o = system(com1)'
                write(10, *) 'h1 = system(com2)'
                write(10, *) 'h2 = system(com3)'
                write(10, *) 'or = 0.74*0.4'
                write(10, *) 'hr = 0.46*0.4'
                write(10, *) 'do for [i=1:tstep] {'
                write(10, *) 'ox = real(word(o, 3*i-2))'
                write(10, *) 'oy = real(word(o, 3*i-1))'
                write(10, *) 'oz = real(word(o, 3*i))'
                write(10, *) 'h1x = real(word(h1, 3*i-2))'
                write(10, *) 'h1y = real(word(h1, 3*i-1))'
                write(10, *) 'h1z = real(word(h1, 3*i))'
                write(10, *) 'h2x = real(word(h2, 3*i-2))'
                write(10, *) 'h2y = real(word(h2, 3*i-1))'
                write(10, *) 'h2z = real(word(h2, 3*i))'
                if ( ilabel == 1 ) then
                        write(10, *) 'set title sprintf("Time = %.2e\n Imaginary frequency", i * dt)'
                else
                        write(10, *) 'set title sprintf("Time = %.2e", i * dt)'
                end if
                write(10, *) 'set urange [0:pi]'
                write(10, *) 'set vrange [0:2*pi]'
                write(10, *) 'splot or*Fx(u,v)+ox, or*Fy(u,v)+oy, or*Fz(u,v)+oz linetype rgbcolor "0xFE0300" notitle, \'
                write(10, *) 'hr*Fx(u,v)+h1x, hr*Fy(u,v)+h1y, hr*Fz(u,v)+h1z linetype rgbcolor "0xFFCCCC" notitle, \'
                write(10, *) 'hr*Fx(u,v)+h2x, hr*Fy(u,v)+h2y, hr*Fz(u,v)+h2z linetype rgbcolor "0xFF7577" notitle'
                write(10, *) '}'

                close(10)

                call system('gnuplot "anime/plot.plt"')

        end subroutine



end program

