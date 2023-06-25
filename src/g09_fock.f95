!**********************************************************************!
!   Read in the fock matrix from a Gaussian output file                !
!**********************************************************************!
subroutine g09_fock( logf )
    use g09_commonvar
    implicit none

    character(*), intent(in) :: logf
    integer, parameter :: fno = 510
    logical exists, read_flag
    character(100) read_line
    integer col(5), row, i, j

    ! check that the log file exists
    inquire( file = logf, exist = exists )
    if ( .not. exists ) then
        print'(3a)', ' Log file file ', logf, ' not found.'
        stop
    end if

    ! open the log file and read the fock matrix
    open( unit = fno, file = logf, status = 'old', &
          action = 'read' )

    ! allocate space for the fock matrix
    if ( allocated ( fock ) ) deallocate( fock )
    allocate( fock( g09_task_numBasisFunctions,          &
                    g09_task_numBasisFunctions ) )

    ! read the fock matrix
    ! we just want the one from the last iteration.
    ! so read in all, but keep overwriting until you 
    ! reach the end of the file
    read_flag = .false.
    do
        read( fno, '(a)', end = 101 ) read_line
        if ( index( read_line , 'Fock matrix (alpha):' ) .ne. 0 ) then
            read_flag = .true.
            do 
                ! first read the column numbers. The formatting in the
                ! log file requires a little extra work to read these
                ! values
                col = 0
                read( fno, '(9x,5(i12,2x))', end = 102 ) col(:)
                102 continue
                ! read the info in each row
                do
                    read ( fno, '(i8, 5D14.6)', end=103 ) row,   &
                        (fock( row, col(i) ), i = 1,           &
                             maxval( col(:) - col(1) ) + 1 )
                    103 continue
                    if ( row == g09_task_numBasisFunctions ) exit
                end do
                if ( maxval( col ) == g09_task_numBasisFunctions ) exit
            end do
        end if
    end do
    101 continue
    if ( read_flag .eqv. .false. ) goto 1001

    ! close the file
    close( fno )

    ! for convience fill in the upper diagonal too
    do i = 1, g09_task_numBasisFunctions
    do j = i+1, g09_task_numBasisFunctions
        fock( i, j ) = fock( j, i )
    end do
    end do

    ! return cleanly
    return

    ! print errors if info not found and abort
    1001 continue
    print'(2a)', ' Fock matrix not found in file ', logf
    print'(2a)', ' Try running Gausion 09 with IOP(5/33=3)'
    stop
  
end subroutine
