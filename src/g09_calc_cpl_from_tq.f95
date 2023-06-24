!**********************************************************************!
!   Calculate the charge transfer integrals from the output of a       !
!   Gaussian calculation.                                              !
!**********************************************************************!
program g09_tq
    use g09_commonvar
    implicit none

    character(200)          :: tqf, fout
    integer, parameter      :: fno = 103
    real(8), allocatable    :: tq(:)
    integer                 :: npoints(3)
    real(8)                 :: dx, startdx(3)
    
    ! initialize the program, read options and input file
    call init( tqf, fout, npoints, startdx, dx)
    
    ! scan phase space and caluclate couplings
    call calc_cpl_from_tq( tqf, fout, npoints, startdx, dx )


end program


!**********************************************************************!
!   Calculate the coupling from two tq files                           ! 
!**********************************************************************!
subroutine calc_cpl_from_tq( tqf, fout, npoints, startdx, dx )
    use g09_commonvar
    implicit none
    
    character(100), intent(in) :: tqf, fout
    integer :: fno = 101
    integer, intent(in) :: npoints(3)
    real(8), intent(in) :: startdx(3), dx
    type( g09_atom_type ), allocatable :: g09_atom1(:), g09_atom2(:)
    integer p1, p2, numAtoms1
    real(8), allocatable :: tq1(:), tq2(:)
    real(8) tqcpl, tdmcpl, mu1(3), mu2(3), com1(3), com2(3), r(3), rmag
    integer x, y, z

    ! read the tq file
    open( unit = fno, file = trim(tqf), status = 'old', action = 'read' )
    read( fno, '(10X, i4)' ) numAtoms1
    allocate( g09_atom1( numAtoms1 ), tq1( numAtoms1 ) )
    read( fno, '(X)' )
    do p1 = 1, numAtoms1
        read( fno, '(5X, 4(f14.7, X))' ) g09_atom1(p1)%x, &
            g09_atom1(p1)%y, g09_atom1(p1)%z, tq1(p1)
    end do
    close( fno )

    ! the second molecule will just be the same as the first
    allocate( g09_atom2( numAtoms1 ), tq2( numAtoms1 ) )
    do p2 = 1, numAtoms1
        tq2(p2) = tq1(p2)
    end do
    close( fno )

    open( unit = fno, file = trim(fout), status = 'new', action = 'write' )
    write( fno, '(2a)' ) "Input file: ", tqf
    write( fno, '(a)' ) "dx (A), dy (A), dz (A), tqcpl (cm-1), tdmcpl (cm-1)"
    ! scan over slipspace
    do x = 1, npoints(1)
    do y = 1, npoints(2)
    do z = 1, npoints(3)
        ! shift the second molecule -- make it based on position of the first
        do p2 = 1, numAtoms1
            g09_atom2(p2)%x = g09_atom1(p2)%x + startdx(1) + (x-1)*dx
            g09_atom2(p2)%y = g09_atom1(p2)%y + startdx(2) + (y-1)*dx
            g09_atom2(p2)%z = g09_atom1(p2)%z + startdx(3) + (z-1)*dx
        end do
        !print*, bohr_to_angstrom*(startdx(1) + (x-1)*dx), bohr_to_angstrom*(startdx(2) + (y-1)*dx), &
        !        bohr_to_angstrom*(startdx(3) + (z-1)*dx)

        ! calculate the coupling (in Hartree)
        tqcpl = 0.d0
        do p1 = 1, numAtoms1
        do p2 = 1, numAtoms1
            tqcpl = tqcpl + tq1(p1)*tq2(p2)/ &
              dsqrt( ( g09_atom1(p1)%x - g09_atom2(p2)%x )**2 + &
                     ( g09_atom1(p1)%y - g09_atom2(p2)%y )**2 + &
                     ( g09_atom1(p1)%z - g09_atom2(p2)%z )**2 )
        end do
        end do
   
        ! convert to wavenumber and print out
        tqcpl = tqcpl * hartree_to_cm
        !print'(a,f14.2)', ' Transition charge coupling (cm-1): ', tqcpl

        ! also calculate according to the point dipole approximation
        ! for comparison

        ! first calculate the transition dipole moments and center of mass
        ! really this isnt the center of mass since the coordinates arent
        ! weighted by mass, but...
        mu1 = 0.d0
        com1 = 0.d0
        do p1 = 1, numAtoms1
            mu1(1)  = mu1(1)  + g09_atom1(p1)%x*tq1(p1)
            mu1(2)  = mu1(2)  + g09_atom1(p1)%y*tq1(p1)
            mu1(3)  = mu1(3)  + g09_atom1(p1)%z*tq1(p1)
            com1(1) = com1(1) + g09_atom1(p1)%x
            com1(2) = com1(2) + g09_atom1(p1)%y
            com1(3) = com1(3) + g09_atom1(p1)%z
        end do
        com1 = com1 / (1.d0 * numAtoms1 )

        mu2 = 0.d0
        com2 = 0.d0
        do p2 = 1, numAtoms1
            mu2(1)  = mu2(1)  + g09_atom2(p2)%x*tq2(p2)
            mu2(2)  = mu2(2)  + g09_atom2(p2)%y*tq2(p2)
            mu2(3)  = mu2(3)  + g09_atom2(p2)%z*tq2(p2)
            com2(1) = com2(1) + g09_atom2(p2)%x
            com2(2) = com2(2) + g09_atom2(p2)%y
            com2(3) = com2(3) + g09_atom2(p2)%z
        end do
        com2 = com2 / (1.d0 * numAtoms1 )
     
        ! the displacement unit vector and the magnitude
        r = com2 - com1
        rmag = dsqrt(sum(r**2))
        r = r/rmag
    
        ! the transition dipole coupling
        tdmcpl = ( mu1(1)*mu2(1) + mu1(2)*mu2(2) + mu1(3)*mu2(3) - &
            3.d0 * ( mu1(1) * r(1) + mu1(2) * r(2) + mu1(3) * r(3) ) * &
                    ( mu2(1) * r(1) + mu2(2) * r(2) + mu2(3) * r(3) ) ) / &
                    rmag**3
        tdmcpl = tdmcpl * hartree_to_cm
        !print'(a,f14.2)', ' Transition dipole coupling (cm-1): ', tdmcpl

        write( fno, '(5(f14.4,","))' ) bohr_to_angstrom*(startdx(1) + (x-1)*dx), &
                               bohr_to_angstrom*(startdx(2) + (y-1)*dx), &
                               bohr_to_angstrom*(startdx(3) + (z-1)*dx), &
                               tqcpl, &
                               tdmcpl 
    end do
    end do
    end do


end subroutine


!**********************************************************************!
!   Initialize the program. Read command line options and input file   ! 
!**********************************************************************!
subroutine init( tqf, fout, npoints, startdx, dx )
    use g09_commonvar
    implicit none

    character(200), intent(out) :: tqf, fout
    real(8), intent(out) :: dx, startdx(3)
    integer, intent(out) :: npoints(3)
    integer nargs, narg, ios, line, pos
    integer, parameter :: fno = 67, fno2 = 68
    character(32) arg, fin, label, fxyz, task
    character(100) buff, emethod
    logical exists, makeinput

    ! check if any command line arguments are found
    nargs = command_argument_count()

    fin=''
    fout=''

    if ( nargs > 0 ) then
        narg = 1
        do while ( narg <= nargs )
            call get_command_argument( narg, arg )
            arg = adjustl(arg)
            select case ( arg )
                case('--help', '-h')
                    call print_help()
                    stop
                case('-i')
                    narg = narg + 1
                    if ( narg > nargs ) then
                        print*, 'No input file given with option -i'
                        call print_help()
                        stop
                    end if
                    ! get the name of the input file
                    call get_command_argument(narg, fin)
                    fin = adjustl(fin)
                    ! check if it exists
                    inquire( file = trim( fin ), exist = exists )
                    if ( .not. exists ) then
                        print'(3a)', ' Input file ', trim(adjustl(fin)), &
                                     ' does not exist.'
                        call print_help()
                        stop
                    end if
                case('-o')
                    narg = narg + 1
                    if ( narg > nargs ) then
                        print*, 'No output file given with option -o'
                        call print_help()
                        stop
                    end if
                    ! get the name of the input file
                    call get_command_argument(narg, fout)
                    fout = adjustl(fout)
                case default
                    print'(3a)', ' Option ', adjustl(arg), ' unknown'
                    call print_help()
                    stop
            end select
            narg = narg + 1
        end do
    else
        call print_help()
        stop
    end if

    if ( fin == '' ) then
        print*, ' No input file was given. Cannot proceed'
        call print_help()
        stop
    else
    ! read the input file
        print* 
        print*, 'Reading input file: ', trim(adjustl(fin))
        open( unit = fno, file = fin, action = 'read' )
        line = 0
        ios = 0
        do while ( ios == 0 )
            ! read the file line by line
            read( fno, '(a)', iostat = ios ) buff
            ! check in out status is ok
            if ( ios == 0 ) then
                ! increase line number
                line = line + 1
                ! find the position of the separator
                pos = scan( buff, ' ' )
                ! store the label
                label = buff( 1:pos )
                ! store the parameter
                buff = buff( pos + 1 : )
                ! ignore comment lines
                if ( label(1:1) ==  '#' ) cycle

                ! set the parameters
                select case (label)
                    case('tqf')
                        read( buff, *, iostat = ios ) tqf
                        print*, 'tqf: ', trim(tqf)
                    case('dx')
                        read( buff, *, iostat = ios ) dx
                        print*, 'dx : ', dx
                    case('npoints')
                        read( buff, *, iostat = ios ) npoints(1), npoints(2), npoints(3)
                        print*, 'npoints(1-3) : ', npoints(1), npoints(2), npoints(3)
                    case('startdx')
                        read( buff, *, iostat = ios ) startdx(1), startdx(2), startdx(3)
                        print*, 'startdx(1-3) : ', startdx(1), startdx(2), startdx(3)
                    case default
                        print*, 'Label ', trim(adjustl(label)), &
                                ' unknown'
                end select
            end if
        end do
        close( fno )
        ! convert distanes to bohr
        dx = dx / bohr_to_angstrom
        startdx(1) = startdx(1) / bohr_to_angstrom
        startdx(2) = startdx(2) / bohr_to_angstrom
        startdx(3) = startdx(3) / bohr_to_angstrom
        print*  
    end if

end subroutine



!**********************************************************************!
!   Print help info for user                                           ! 
!**********************************************************************!
subroutine print_help()
    implicit none
    
    print*  
    print*, ' g09_calc_cpl_from_tq version 1.0'
    print*  
    print*, ' This  program  calculates  the coulombic coupling from     '
    print*, ' transition charges calculated from g09_tq assuming two     '
    print*, ' identical molecules experincing some displacement.         '
    print* 
    print*, ' Usage:'
    print* 
    print*, ' -- help, -h: print this message'
    print*, ' -i file    : specifies the input file containing the  name '
    print*, '              of the tq file from g09_tq and the parameters '
    print*, '              for the scan.                                 '
    print*, ' -o file    : specifies the output file'
    print*  
    print* 

end subroutine
