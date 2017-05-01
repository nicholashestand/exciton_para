!**********************************************************************!
!   Calculate the charge transfer integrals from the output of a       !
!   Gaussian calculation.                                              !
!**********************************************************************!
program g09_tq
    use g09_commonvar
    implicit none

    character(200) :: fch, logf, fout
    integer, parameter :: fno = 103
    real(8), allocatable :: tq(:)
    integer p, b, c, i, j, estate
    integer, allocatable :: basis_to_atom_map(:)
    logical warn
    real(8) tq_mu(3), tq_sum


    ! initialize the program, read options and input file
    call tq_init(fch, logf, fout, estate )
    
    ! read in info from gaussian files
    print*, '>> Getting info about the molecule'
    call g09_calcInfo( fch )
    call g09_mocoeff( fch )
    call g09_overlap( logf )
    call g09_cicoeff( logf, estate )

    ! map the ao basis functions to the atoms
    allocate( basis_to_atom_map(g09_task_numBasisFunctions) )
    c = 0
    do p = 1, g09_task_numAtoms
    do b = 1, g09_atom(p)%basisNum
        c = c + 1
        basis_to_atom_map(c) = p
    end do
    end do

    ! calculate the transition charges
    allocate(tq(g09_task_numAtoms))
    tq = 0.d0

    ! loop over all basis functions
    warn = .true.
    do b = 1, g09_task_numBasisFunctions
        p = basis_to_atom_map(b)
    do c = 1, g09_task_numBasisFunctions

        ! truncate if the overlap is zero
        if ( overlap(b,c) == 0.d0 ) cycle

        ! loop over all molecular orbitals
        do i = 1, g09_task_numBasisFunctions
        do j = 1, g09_task_numBasisFunctions

            ! truncate if the CI coefficients are zero
            if ( cic(i,j) == 0.d0 ) cycle
        
            if ( i < j ) then
                tq(p) = tq(p) + cic(i,j)*moc(b,i)*moc(c,j)*overlap(b,c)
            else if ( warn ) then
                warn = .false.
                print*, '---------------------------------------------'
                print*, 'Warning... cannot handle deexcitations  <- . '
                print*, 'The current CI coefficient is being ignored. '
                print*, 'Try using CIS instead of TD for the  excited '
                print*, 'state calculation.'
                print*, '---------------------------------------------'
            end if
        
        end do
        end do
    end do
    end do

    ! calcualte the tansition dipole moment and sum of transition charges
    tq_mu = 0.d0
    tq_sum = 0.d0
    do p = 1, g09_task_numAtoms
        tq_sum = tq_sum + tq(p) 
        tq_mu(1) = tq_mu(1) + g09_atom(p)%x*tq(p)
        tq_mu(2) = tq_mu(2) + g09_atom(p)%y*tq(p)
        tq_mu(3) = tq_mu(3) + g09_atom(p)%z*tq(p)
    end do

    open( unit = fno, file = trim(fout), status = 'new', action = 'write')
    write( fno, * ) 'Atomic Number, X (bohr), Y(bohr), Z(bohr), TQ (au)'
    do p = 1, g09_task_numAtoms 
        write( fno, '(i4,",",4(f14.7,","))' ) g09_atom(p)%atomicNum, &
            g09_atom(p)%x, g09_atom(p)%y, g09_atom(p)%z, tq(p)
    end do
    write( fno, * ) 'tq sum,', tq_sum
    write( fno, '(a,",",3(f14.7,","))' ) ' tq dipole (au)', tq_mu(1), &
                                         tq_mu(2), tq_mu(3)

    tq_mu = tq_mu * au_to_debye
    write( fno, '(a,",",3(f14.7,","))' ) ' tq dipole (debye)', tq_mu(1), &
                                         tq_mu(2), tq_mu(3)

    print*, '>> Wrote output to file ', trim(adjustl(fout))
    print*, 

end program



!**********************************************************************!
!   Initialize the program. Read command line options and input file   ! 
!**********************************************************************!
subroutine tq_init(fch, logf, fout, estate )
    implicit none

    character(200), intent(out) :: fch, logf, fout
    integer, intent(out) :: estate
    integer nargs, narg, ios, line, pos
    integer, parameter :: fno = 67, fno2 = 68
    character(32) arg, fin, label, fxyz, task
    character(100) buff, emethod
    logical exists, makeinput

    makeinput = .false.
    fin = ''
    fout = ''
    fxyz = ''
   

    ! check if any command line arguments are found
    nargs = command_argument_count()

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
                case('-m')
                    makeinput = .true. 
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
        print*,
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
                    case('emethod')
                        read( buff, *, iostat = ios ) emethod
                        print*, 'emethod : ', trim(adjustl(emethod))
                    case('task')
                        read( buff, *, iostat = ios ) task
                        print*, 'task : ', trim(adjustl(task))
                    case('fch')
                        read( buff, *, iostat = ios ) fch
                        print*, 'fch: ', trim(adjustl(fch))
                    case('log')
                        read( buff, *, iostat = ios ) logf
                        print*, 'log: ', trim(adjustl(logf))
                    case('xyz_file')
                        read( buff, *, iostat = ios ) fxyz
                        print*, 'fxyz: ', trim(adjustl(fxyz))
                    case('estate')
                        read( buff, *, iostat = ios ) estate
                        print*, 'estate: ', estate
                    case default
                        print*, 'Label ', trim(adjustl(label)), &
                                ' unknown'
                end select
            end if
        end do
        close( fno )
        print*, 
    end if

    if ( makeinput ) then
        ! set default method and basis if they are not set in the input file
        if ( emethod == '' ) emethod = 'CIS(Nstates=1,Root=1)/cc-pVTZ'
        
        ! make sure the xyz files exists
        inquire( file = trim(fxyz), exist = exists )
        if ( .not. exists ) then
            print*, 'The file ', trim(adjustl(fxyz)), &
                    ' does not exist. Aborting...'
            stop
        end if
 
        ! setup the gaussian calculation
        open( unit = fno, file = trim(adjustl(task))//'.gjf', &
              action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))//'.chk'
        write( fno, * ) '%NProcShared=4'
        write( fno, * ) '%mem=1GB'
        write( fno, * ) '# '//trim(emethod)//' SP IOP(9/40=5) IOP(3/33=1)' &
                        //' Pop=Full Symmetry=None'
        write( fno, * )
        write( fno, * ) trim(adjustl(task))
        write( fno, * )
        write( fno, * ) '0 1'
        ! write the coordinates from the xyz file
        open( unit = fno2, file = trim(fxyz), action = 'read' )
        do
            read( fno2, '(a)', end=105) buff
            write( fno, * ) buff
        end do
105     continue
        close( fno2 )
        close( fno )
    
        print*, 'Successfuly made the G09 input files'
        stop
    end if

end subroutine



!**********************************************************************!
!   Print help info for user                                           ! 
!**********************************************************************!
subroutine print_help()
    implicit none

    print*, 
    print*, ' g09_tq version 1.0'
    print*, 
    print*, ' This  program  calculates  the transition  charges  for  a '
    print*, ' specified excited state  of a  molecule.  The  calculation '
    print*, ' uses output from a Gaussian 09 calculation.'
    print*,
    print*, ' Usage:'
    print*,
    print*, ' -- help, -h: print this message'
    print*, ' -i file    : specifies the input file containing the names '
    print*, '              of the formatted checkpoint and log files. If '
    print*, '              the -m option is also given, the  input  file '
    print*, '              is read for information about the Gaussian 09 '
    print*, '              calculation to add to the input files.'
    print*, ' -o file    : specifies the output file'
    print*, ' -m         : makes one Gaussian 09 input files from  a xyz '
    print*, '              files of atomic coordinates.'
    print*, 
    print*, ' To  use  this  program to calculate the transition charges '
    print*, ' you  must  first  run  an  excited  state  calculation  in '
    print*, ' Gaussian 09 of the molecule of interest. The -m  option of '
    print*, ' this program can  be  used  to  generate  the  appropriate '
    print*, ' Gaussian 09 input files if xyz  coordinate files  for  the '
    print*, ' molecule of interest is provided.'
    print*, 
    print*, ' The  resulting  checkpoint  file   must  be  converted  to '
    print*, ' formatted checkpoint files using the Gaussian 09 formcheck '
    print*, ' utility before  it  can be read by this program.'
    print*,

end subroutine


