!**********************************************************************!
!   Get information about a Gaussian 09 calculation from the           !
!   formatted checkpoint file. If not already, all quantities are      !
!   converted to atomic units                                          !
!**********************************************************************!
subroutine g09_calcInfo( fch )
    use g09_commonvar
    implicit none

    character(*), intent(in) :: fch
    integer, parameter :: fno = 509 
    logical exists
    character(100) read_line
    integer i

    ! check that the formatted checkpoint file exists
    inquire( file = fch, exist = exists )
    if ( .not. exists ) then
        print'(3a)', ' Formatted checkpoint file ', trim(adjustl(fch)), &
                     ' not found.'
        stop
    end if

    ! open the formatted checkpoint file and read calculation info
    open( unit = fno, file = fch, status = 'old', &
          action = 'read' )
    read( fno, '(a)' ) g09_task_title
    read( fno, '(a10, a30, 30x, a30)' ) g09_task_type, g09_task_method, g09_task_basis
    read( fno, '(49x, i12)') g09_task_numAtoms
    ! read the number of electrons
    do 
        read( fno, '(a)', end = 1002 ) read_line
        if ( index( read_line, 'Number of electrons' ) .ne. 0 ) then
            read( read_line, '(49x, i12)' ) g09_task_numElectrons
            exit
        end if
    end do
    ! read the number of basis functions
    do 
        read( fno, '(a)', end = 1003 ) read_line
        if ( index( read_line, 'Number of basis functions' ) .ne. 0 ) then
            read( read_line, '(49x, i12)' ) g09_task_numBasisFunctions
            exit
        end if
    end do
    ! read the number of independent basis functions
    do 
        read( fno, '(a)', end = 1011 ) read_line
        if ( index( read_line, 'Number of independent functions' ) .ne. 0 ) then
            read( read_line, '(49x, i12)' ) g09_task_numBasisFunctionsUsed
            exit
        end if
    end do

    ! allocate space for the atom information and read it in
    if ( allocated( g09_atom )) deallocate(g09_atom)
    allocate( g09_atom( g09_task_numAtoms ) )

    ! read the atomic numbers
    do
        read( fno, '(a)', end = 1004 ) read_line
        if ( index( read_line , 'Atomic numbers' ) .ne. 0 ) then
            read( fno, '(6i12)' ) g09_atom(:)%atomicNum
            exit
        end if
    end do

    ! read the cartesian coordinates
    do
        read( fno, '(a)', end = 1005 ) read_line
        if ( index( read_line , 'Current cartesian coordinates' ) .ne. 0 ) then
            read( fno, '(5E16.8)' ) (g09_atom(i)%x, g09_atom(i)%y,          &
                                     g09_atom(i)%z, i=1, g09_task_numAtoms )
            exit
        end if
    end do

    ! read the atomic masses
    do
        read( fno, '(a)', end = 1010 ) read_line
        if ( index( read_line , 'Real atomic weights' ) .ne. 0 ) then
            read( fno, '(5E16.8)' ) g09_atom(:)%mass
            exit
        end if
    end do
    ! convert from amu to electron mass
    g09_atom(:)%mass = g09_atom(:)%mass * AMU_to_electronMass



    ! calculate the number of degrees of freedom
    g09_task_dof = 3*g09_task_numAtoms



    ! allocate space for and read the atomic orbital basis information
    ! read shell types
    do 
        read( fno, '(a)', end = 1006 ) read_line
        if ( index( read_line, 'Shell types' ) .ne. 0 ) then
            read( read_line, '(49x, i12)' ) g09_task_numAOShells
            if ( allocated( g09_AOShell ) ) deallocate( g09_AOShell )
            allocate( g09_AOShell( g09_task_numAOShells ) )
            read( fno, '(6i12)' ) g09_AOShell(:)%typ
            exit
        end if
    end do
    ! read shell to atom map
    do 
        read( fno, '(a)', end = 1007 ) read_line
        if ( index( read_line, 'Shell to atom map' ) .ne. 0 ) then
            read( fno, '(6i12)' ) g09_AOShell(:)%map
            exit
        end if
    end do
    ! assign number of basis functions to each atom
    g09_atom%basisNum = 0
    do i = 1, g09_task_numAOShells
        g09_atom( g09_AOShell(i)%map )%basisNum =       &
             g09_atom( g09_AOShell(i)%map )%basisNum    &
            +g09_AOShellDef(g09_AOShell(i)%typ)
    end do

    ! read the total energy
    do 
        read( fno, '(a)', end = 1008 ) read_line
        if ( index( read_line, 'Total Energy' ) .ne. 0 ) then
            read( read_line, '(49x, E22.15)' ) g09_task_totalEnergy
            exit
        end if
    end do
 


    ! close the formatted checkpoint file   
    close( fno )

    ! print out calculation info
    print'(a)',         ' Gaussian calculation details: '
    print'(a)',         '    Input File:            '//trim(g09_task_title)
    print'(a)',         '    Method:                '//trim(g09_task_method)
    print'(a)',         '    Basis:                 '//trim(g09_task_basis)
    print'(a,i12)',     '    Atoms:                 ', g09_task_numAtoms
    print'(a,i12)',     '    Electrons:             ', g09_task_numElectrons
    print'(a,i12)',     '    Basis Functions:       ', g09_task_numBasisFunctions
    print'(a,i12)',     '    Basis Functions Used:  ', g09_task_numBasisFunctionsUsed
    print'(a,E12.6)',   '    Energy (hartree):      ', g09_task_totalEnergy

    ! return cleanly
    return

    ! print errors if info not found and abort
    1002 continue
    print'(2a)', 'Number of electrons not found in file ', fch
    stop
    1003 continue
    print'(2a)', 'Number of basis functions not found in file ', fch
    stop
    1004 continue
    print'(2a)', 'Atomic numbers not found in file ', fch
    stop
    1005 continue
    print'(2a)', 'Cartesian coordinates not found in file ', fch
    stop
    1006 continue
    print'(2a)', 'Shell types not found in file ', fch
    stop
    1007 continue
    print'(2a)', 'Shell to atom map not found in file ', fch
    stop
    1008 continue
    print'(2a)', 'Total energy not found in file ', fch
    stop
    1010 continue
    print'(2a)', 'Real atomic weights not found in file ', fch
    stop
    1011 continue
    print'(2a)', 'Number of independent basis functions not found in file ', fch
    stop

end subroutine
