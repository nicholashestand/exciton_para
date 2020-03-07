!**********************************************************************!
!   Read in the moleuclar orbital coefficients from the formatted      !
!   checkpoint file. Remember, these are in a basis of atomic          !
!   orbitals.                                                          !
!**********************************************************************!
subroutine g09_mocoeff( fch )
    use g09_commonvar
    implicit none

    character(*), intent(in) :: fch
    integer, parameter :: fno = 509 
    logical exists
    character(100) read_line

    ! check that the formatted checkpoint file exists
    inquire( file = fch, exist = exists )
    if ( .not. exists ) then
        print'(3a)', ' Formatted checkpoint file ', fch, &
                     ' not found.'
        stop
    end if

    ! open the formatted checkpoint file and read the molecular orbital
    ! coefficient matrix and the molecular orbital energies
    open( unit = fno, file = fch, status = 'old', &
          action = 'read' )

    ! allocate space for the molecular orbital energies
    if ( allocated( moe ) ) deallocate(moe)
    allocate( moe( g09_task_numBasisFunctionsUsed ) )
    do
        read( fno, '(a)', end = 1001 ) read_line
        if ( index( read_line , 'Alpha Orbital Energies' ) .ne. 0 ) then
            read( fno, '(5E16.8)' ) moe(:)
            exit
        end if
    end do

    ! allocate space for the molecular orbital coefficients
    if ( allocated( moc ) ) deallocate( moc )
    allocate( moc( g09_task_numBasisFunctions,g09_task_numBasisFunctionsUsed ) )

    ! read the molecular orbital coefficients
    do
        read( fno, '(a)', end = 1002 ) read_line
        if ( index( read_line , 'Alpha MO coefficients' ) .ne. 0 ) then
            read( fno, '(5E16.8)' ) moc(:,:)
            exit
        end if
    end do

    ! close the file
    close( fno )

    ! return cleanly
    return

    ! print errors if info not found and abort
    1001 continue
    print'(2a)', 'Alpha Orbital Energies not found in file ', fch
    stop
    1002 continue
    print'(2a)', 'Alpha MO coefficients not found in file ', fch
    stop
  
end subroutine
