!**********************************************************************!
!   Calculate the vibronic coupling parameter from the output of a     !
!   Gaussian calculation.                                              !
!**********************************************************************!
program g09_lambda
    use g09_commonvar
    implicit none

    character(200) :: fch_g, fch_e,         &
                      fch_a, fch_c, fout
    type (g09_atom_type), allocatable :: g09_atom_g(:),     &
                                         g09_atom_e(:),     &
                                         g09_atom_a(:),     &
                                         g09_atom_c(:)
    real(8), allocatable :: freq(:), normalModes(:,:),      &
                            lambda_ge(:), lambda_ga(:),     &
                            lambda_gc(:)
    real(8) :: lbnd, ubnd, hreff, freqeff
    integer i
    integer, parameter :: fno = 607


    ! initialize the program, read options and input file
    call lambda_init(fch_g, fch_e, fch_a, fch_c, fout, lbnd, ubnd)
    
    ! read in the data from the Guassian files
    print*, '>> Getting info about the ground (reference) state.'
    call g09_calcInfo(fch_g) ! returns several variables from g09_commonvar
    call g09_hessian(fch_g)  ! returns the hessian from the ground state
    allocate( g09_atom_g(g09_task_numAtoms),            &
              g09_atom_e(g09_task_numAtoms),            &
              g09_atom_a(g09_task_numAtoms),            &
              g09_atom_c(g09_task_numAtoms),            &
              freq(g09_task_dof),                       &
              normalmodes(g09_task_dof, g09_task_dof))
    g09_atom_g = g09_atom    ! store the coordinates associated with 
                             ! the ground state


    ! calculate the vibrational frequencies
    call calc_vib_freq(g09_atom_g, normalmodes, freq)

    ! calculate the vibronic coupling parameters
    allocate( lambda_ge(g09_task_dof),      &
              lambda_ga(g09_task_dof),      &
              lambda_gc(g09_task_dof) )
    lambda_ge = -1.d10
    lambda_ga = -1.d10
    lambda_gc = -1.d10

    ! ground to excited
    if ( fch_e .ne. '' ) then
        print*, '>> Getting info about the excited (deformed) state.'
        ! get info about the excited state calculation
        call g09_calcInfo(fch_e)
        g09_atom_e = g09_atom
        ! calculate vibronic coupling parameters
        call calc_VibronicCoupling( g09_atom_g, g09_atom_e,  &
                                    normalmodes, freq,       &
                                    lambda_ge )
    end if

    ! ground to anion
    if ( fch_a .ne. '' ) then
        print*, '>> Getting info about the anion (deformed) state.'
        ! get info about the anion calculation
        call g09_calcInfo(fch_a)
        g09_atom_a = g09_atom
        ! calculate vibronic coupling parameters
        call calc_VibronicCoupling( g09_atom_g, g09_atom_a,  &
                                    normalmodes, freq,       &
                                    lambda_ga )
    end if

    ! ground to cation
    if ( fch_c .ne. '' ) then
        print*, '>> Getting info about the cation (deformed) state.'
        ! get info about the cation calculation
        call g09_calcInfo(fch_c)
        g09_atom_c = g09_atom
        ! calculate vibronic coupling parameters
        call calc_VibronicCoupling( g09_atom_g, g09_atom_c,  &
                                    normalmodes, freq,       &
                                    lambda_gc )
    end if

    ! write output to file
    if ( .not. fout == '' ) then
        ! convert frequencies to cm-1
        freq = freq * Hartree_to_cm
        open ( unit = fno, file = fout, status = 'new', action = 'write' )
        write( fno, * ) fout
        write( fno, * ) g09_task_method
        write( fno, * ) g09_task_basis
        if ( fch_e .ne. '' ) then
            write( fno, * )
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            write( fno, '(3a20)' ) 'Mode', 'Frequency (cm-1)', 'lambda'
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            do i = 1, g09_task_dof
                write( fno, '(i20, 2f20.10)' ) i, freq(i), lambda_ge(i)
            end do
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            write( fno, '(a,2f10.1)' ) 'Effective HR range (cm-1):', &
                                        lbnd, ubnd
            hreff = 0.d0
            freqeff = 0.d0
            do i = 1, g09_task_dof
                if ( freq(i) > lbnd .and. freq(i) < ubnd ) then
                    hreff = hreff + lambda_ge(i)**2
                    freqeff = freqeff + freq(i)*lambda_ge(i)**2
                end if
            end do
            write( fno, '(a,f10.2)' ) 'Effective HR Factor: ', hreff
            write( fno, '(a,f10.2)' ) 'Effective HR Frequency: ', &
                                       freqeff/hreff
        end if
        if ( fch_a .ne. '' ) then
            write( fno, * )
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            write( fno, '(3a20)' ) 'Mode', 'Frequency (cm-1)', 'lambda-'
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            do i = 1, g09_task_dof
                write( fno, '(i20, 2f20.10)' ) i, freq(i), lambda_ga(i)
            end do
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            write( fno, '(a,2f10.1)' ) 'Effective HR- range (cm-1):', &
                                        lbnd, ubnd
            hreff = 0.d0
            freqeff = 0.d0
            do i = 1, g09_task_dof
                if ( freq(i) > lbnd .and. freq(i) < ubnd ) then
                    hreff = hreff + lambda_ga(i)**2
                    freqeff = freqeff + freq(i)*lambda_ga(i)**2
                end if
            end do
            write( fno, '(a,f10.2)' ) 'Effective HR- Factor: ', hreff
            write( fno, '(a,f10.2)' ) 'Effective HR- Frequency:', &
                                       freqeff/hreff
        end if
        if ( fch_c .ne. '' ) then
            write( fno, * )
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            write( fno, '(3a20)' ) 'Mode', 'Frequency (cm-1)', 'lambda+'
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            do i = 1, g09_task_dof
                write( fno, '(i20, 2f20.10)' ) i, freq(i), lambda_gc(i)
            end do
            write( fno, '(a60)' ) '***********************************'//&
                                  '*************************'
            write( fno, '(a,2f10.1)' ) 'Effective HR+ range (cm-1):', &
                                        lbnd, ubnd
            hreff = 0.d0
            freqeff = 0.d0
            do i = 1, g09_task_dof
                if ( freq(i) > lbnd .and. freq(i) < ubnd ) then
                    hreff = hreff + lambda_gc(i)**2
                    freqeff = freqeff + freq(i)*lambda_gc(i)**2
                end if
            end do
            write( fno, '(a,f10.2)' ) 'Effective HR+ Factor: ', hreff
            write( fno, '(a,f10.2)' ) 'Effective HR+ Frequency: ', &
                                      freqeff/hreff
        end if
        close( fno ) 
        print*, 
        print*, '>> Wrote output to file: ', trim(adjustl(fout))
        print*,
    end if

    
end program



!**********************************************************************!
!   Calculate the vibrational frequencies from the Hessian             !
!**********************************************************************!
subroutine calc_vib_freq(g09_atom_g, normalModes, freq)
    use g09_commonvar
    implicit none

    type( g09_atom_type ), intent(in) :: g09_atom_g(g09_task_numAtoms)
    real(8), intent(out) :: freq(g09_task_dof)
    real(8), intent(out) :: normalModes( g09_task_dof, g09_task_dof )
    integer i, j, atom1, atom2
    real(8) hessian_mwc(g09_task_dof, g09_task_dof),    &
            rmass(g09_task_dof), kforce(g09_task_dof),  &
            mass

    ! calculate the vibrational frequencies
    print*, '>> Calculating the vibrational frequencies...'
    print*, '>> There are ', g09_task_dof, ' degrees of freedom.'

    ! calculate the mass weighted hessian
    do i = 1, g09_task_dof
    do j = 1, g09_task_dof
        atom1 = ceiling(i/3.d0)
        atom2 = ceiling(j/3.d0)
        mass = dsqrt( g09_atom_g(atom1)%mass * g09_atom_g(atom2)%mass )
        hessian_mwc(i,j) = hessian(i,j)/mass
    end do
    end do

    ! diagonalize the mass weighted hessian to find the vibrational 
    ! normal modes
    call diagonalize( hessian_mwc, g09_task_dof, freq )

    ! save the normal modes to send back
    normalModes = hessian_mwc
    
    ! the frequencies are the square root of the eienvalues of the Hessian,
    ! imaginary frequencies are just made negative, since there are no
    ! imaginary frequencies
    do i = 1, g09_task_dof
        if ( freq(i) < 0 ) then
            freq(i) = -dsqrt(abs(freq(i)))
        else
            freq(i) = dsqrt(freq(i))
        end if
    end do

    ! calculate the reduced mass of each normal mode from the eigenvectors 
    ! of the mass weighted hessian. From here on, I calculate these just
    ! so they can be compared with the results from the Gaussian calculation
    rmass = 0.d0
    do i = 1, g09_task_dof
    do j = 1, g09_task_dof
        atom1 = ceiling(i/3.d0)
        rmass(j) = rmass(j) + hessian_mwc(i,j)**2/g09_atom_g(atom1)%mass
    end do
    end do
    rmass = 1.d0/rmass

    ! calculate the force constant in Hartree**2*massElectron
    kforce = freq**2 * rmass

    ! print out useful info to terminal to compare with
    ! output from gaussian calculation
    print*, '************************************'// &
            '************************************'
    print*, '    Properties of the lowest normal modes'
    print*, '************************************'// &
            '************************************'
    print'(3a20)', 'Frequency', 'Reduced Mass', 'Force Constant'
    print'(3a20)', '(cm-1)', '(amu)', '(mdyne/A)'
    do i = 1, min(20, g09_task_dof)
        print'(3(f20.4))', freq(i) *Hartree_to_cm,            &
                           rmass(i)/AMU_to_electronMass,      &
                           kforce(i) * au_to_mdyneA
    end do
    print*, '************************************'// &
            '************************************'
    
end subroutine



!**********************************************************************!
!   Call the lapack diagonalization function dsyev                     !
!**********************************************************************!
subroutine diagonalize( A, N, W )
    implicit none

    integer, intent(in) :: N
    real(8), intent(inout) :: A(N,N)
    real(8), intent(out) :: W(N)
    character(4) :: jobz = 'v', uplo = 'u'
    real(8), allocatable :: work(:)
    integer lwork, nb, info
    integer, external :: ilaenv

    ! find optimal work space
    nb = ilaenv( 1, 'dsytrd', 'u', N, -1, -1, -1 )
    lwork = (nb+2) * N
    allocate (work(lwork))

    ! ask lapack to diagonalize the matrix
    call dsyev( jobz, uplo, N, A, N, W, work, lwork, info )

end subroutine



!**********************************************************************!
!   Calculate the vibronic coupling parameter                          !
!**********************************************************************!
subroutine calc_VibronicCoupling( g09_atom_ref, g09_atom_displaced, &
                                  normalModes, freq, lambda )
    use g09_commonvar
    implicit none

    type(g09_atom_type), intent(in) ::                 &
                g09_atom_ref( g09_task_numAtoms ),     &
                g09_atom_displaced( g09_task_numAtoms )
    real(8), intent(in) :: freq(g09_task_dof)
    real(8), intent(in) :: normalModes(g09_task_dof, g09_task_dof)
    real(8), intent(out) :: lambda(g09_task_dof)
    real(8) :: displacement(g09_task_dof), Ro(g09_task_dof)
    integer :: nx, atom

    ! find the geometry displacement in mass weighted coordinates
    ! units is bohr sqrt(electronMass)
    do atom = 1, g09_task_numAtoms
        nx = 1 + 3*(atom-1) ! x
        displacement( nx ) = ( g09_atom_displaced( atom )%x - &
                               g09_atom_ref( atom )%x) *      &
                               dsqrt(g09_atom_ref( atom )%mass )
        nx = 2 + 3*(atom-1) ! y
        displacement( nx ) = ( g09_atom_displaced( atom )%y - &
                               g09_atom_ref( atom )%y) *      &
                               dsqrt(g09_atom_ref( atom )%mass )
        nx = 3 + 3*(atom-1) ! z
        displacement( nx ) = ( g09_atom_displaced( atom )%z - &
                               g09_atom_ref( atom )%z) *      &
                               dsqrt(g09_atom_ref( atom )%mass )
    end do

    ! project the displacements onto the normal modes. 
    Ro = matmul( transpose( normalModes ), displacement ) 

    ! the lambda (lambda**2 is the Huang-Rhys factor)
    ! lambda is dimensionless
    lambda = dsqrt(dabs(freq)/2.d0)*Ro 

end subroutine



!**********************************************************************!
!   Initialize the program. Read command line options and input file   ! 
!**********************************************************************!
subroutine lambda_init(fch_g, fch_e, fch_a, fch_c, fout, lbnd,  ubnd)
    implicit none

    character(200), intent(out) :: fch_g, fch_e, &
                                   fch_a, fch_c, fout
    real(8), intent(out) :: lbnd, ubnd
    integer nargs, narg, ios, line, pos
    integer, parameter :: fno = 67, fno2 = 68
    character(32) arg, fin, label, fxyz, task
    character(100) buff
    logical exists, makeinput
    character(64) method, emethod

    makeinput = .false.
    fin = ''
    fout = ''

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
                    case('method')
                        read( buff, *, iostat = ios ) method
                        print*, 'method : ', trim(adjustl(method))
                    case('emethod')
                        read( buff, *, iostat = ios ) emethod
                        print*, 'emethod : ', trim(adjustl(emethod))
                    case('task')
                        read( buff, *, iostat = ios ) task
                        print*, 'task : ', trim(adjustl(task))
                    case('ground_fch')
                        read( buff, *, iostat = ios ) fch_g
                        print*, 'fch_g: ', trim(adjustl(fch_g))
                    case('excited_fch')
                        read( buff, *, iostat = ios ) fch_e
                        print*, 'fch_e: ', trim(adjustl(fch_e))
                    case('anion_fch')
                        read( buff, *, iostat = ios ) fch_a
                        print*, 'fch_a: ', trim(adjustl(fch_a))
                    case('cation_fch')
                        read( buff, *, iostat = ios ) fch_c
                        print*, 'fch_c: ', trim(adjustl(fch_c))
                    case('xyz_file')
                        read( buff, *, iostat = ios ) fxyz
                        print*, 'fxyz: ', trim(adjustl(fxyz))
                    case('lower_bound')
                        read( buff, *, iostat = ios ) lbnd
                        print*, 'lbnd: ', lbnd
                    case('upper_bound')
                        read( buff, *, iostat = ios ) ubnd
                        print*, 'ubnd: ', ubnd
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
        ! set default method and basis if they are not set 
        ! in the input file
        if ( method == '' ) method = 'HF/6-31G'
        if ( emethod == '' ) emethod = 'CIS(nstates=1,root=1)/6-31G'
        
        ! make sure the xyz file exists
        inquire( file = trim(fxyz), exist = exists )
        if ( .not. exists ) then
            print*, 'The file ', trim(adjustl(fxyz)),&
                    ' does not exist. Aborting...'
            stop
        end if
 
        ! the ground state calculation
        open( unit = fno, file = trim(adjustl(task))//'_g.gjf', &
                          action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))//'_g.chk'
        write( fno, * ) '%NProcShared=4'
        write( fno, * ) '# opt=tight freq=(savenormalmodes,hpmodes) '//&
                         trim(method)//&
                         ' Symmetry=Loose int=ultrafine MaxDisk=2GB'
        write( fno, * )
        write( fno, * ) trim(adjustl(task))//'_g'
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

        ! the excited state calculation
        open( unit = fno, file = trim(adjustl(task))//'_e.gjf', &
                          action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))//'_e.chk'
        write( fno, * ) '%NProcShared=4'
        write( fno, * ) '# opt=tight '//trim(emethod)//&
                        ' Symmetry=Loose int=ultrafine MaxDisk=2GB'
        write( fno, * )
        write( fno, * ) trim(adjustl(task))//'_e'
        write( fno, * )
        write( fno, * ) '0 1'
        ! write the coordinates from the xyz file
        open( unit = fno2, file = trim(fxyz), action = 'read' )
        do
            read( fno2, '(a)', end=106) buff
            write( fno, * ) buff
        end do
106     continue
        close( fno2 )
        close( fno )
        
        ! the anion state calculation
        open( unit = fno, file = trim(adjustl(task))//'_a.gjf',&
                          action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))//'_a.chk'
        write( fno, * ) '%NProcShared=4'
        write( fno, * ) '# opt=tight '//trim(method)//&
                        ' Symmetry=Loose int=ultrafine MaxDisk=2GB'
        write( fno, * )
        write( fno, * ) trim(adjustl(task))//'_a'
        write( fno, * )
        write( fno, * ) '-1 2'
        ! write the coordinates from the xyz file
        open( unit = fno2, file = trim(fxyz), action = 'read' )
        do
            read( fno2, '(a)', end=107) buff
            write( fno, * ) buff
        end do
107     continue
        close( fno2 )
        close( fno )

        ! the cation state calculation
        open( unit = fno, file = trim(adjustl(task))//'_c.gjf', &
                          action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))//'_c.chk'
        write( fno, * ) '%NProcShared=4'
        write( fno, * ) '# opt=tight '//trim(method)//&
                        ' Symmetry=Loose int=ultrafine MaxDisk=2GB'
        write( fno, * )
        write( fno, * ) trim(adjustl(task))//'_c'
        write( fno, * )
        write( fno, * ) '1 2'
        ! write the coordinates from the xyz file
        open( unit = fno2, file = trim(fxyz), action = 'read' )
        do
            read( fno2, '(a)', end=108) buff
            write( fno, * ) buff
        end do
108     continue
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
    print*, ' g09_lambda version 1.0'
    print*, 
    print*, ' This program calculates the vibronic coupling parameter for '
    print*, ' each  vibrational mode of a molecule  using  the  formatted '
    print*, ' checckpoint file from Gaussian09 output.'
    print*,
    print*, ' Usage:'
    print*,
    print*, ' -- help, -h: print this message'
    print*, ' -i file    : specifies the input file containing the names '
    print*, '              of the formatted checkpoing files. If the  -m '
    print*, '              option is also given, it reads the input file '
    print*, '              for information about g09 calculation to  add '
    print*, '              to the input file.'
    print*, ' -o file    : specifies the output file.'
    print*, ' -m         : makes four gaussian input  files from  an xyz '
    print*, '              file of atomic coordinates.  The  four  files '
    print*, '              are for the ground, excited, anion and cation '
    print*, '              optimizations. The options are  set  so  that '
    print*, '              the resulting fch files can  be used directly '
    print*, '              with this program.  Note that  the  formatted '
    print*, '              checkpoint  file  needs  to  be  made  by the '
    print*, '              Gaussian formchk utility.'
    print*, 
    print*, ' To use this  program to  calculate the  vibronic  coupling '
    print*, ' parameters you must  first  run  a  geometry  optimization '
    print*, ' in  Gaussian 09 of the  reference  (ground)  and displaced '
    print*, ' (excited) states. The reference state must also  include a '
    print*, ' frequency calculation so that the Hessian gets  written to '
    print*, ' the checkpoint file. The -m option of this program can  be '
    print*, ' used to generate the appropriate Gaussian 09  input  files '
    print*, ' if a xyz coordinate file for the molecule of  interest  is '
    print*, ' provided.'
    print*, 
    print*, ' The  resulting  checkpoing  files  must  be  converted  to '
    print*, ' formatted checkpoint files using the Gaussian 09 formcheck '
    print*, ' utility before they can be read by this program.'
    print*,

end subroutine


