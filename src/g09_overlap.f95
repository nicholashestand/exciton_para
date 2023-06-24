!**********************************************************************!
!   Read in the overlap matrix from a Gaussian output file             !
!**********************************************************************!
subroutine g09_overlap( logf )
    use g09_commonvar
    implicit none

    character(*), intent(in) :: logf
    integer, parameter :: fno = 509 
    logical exists
    character(100) read_line
    integer col(5), row, i, j

    ! check that the log file exists
    inquire( file = logf, exist = exists )
    if ( .not. exists ) then
        print'(3a)', ' Log file file ', logf, ' not found.'
        stop
    end if

    ! open the log file and read the overlap matrix
    open( unit = fno, file = logf, status = 'old', &
          action = 'read' )

    ! allocate space for the overlap matrix
    if ( allocated ( overlap ) ) deallocate( overlap )
    allocate( overlap( g09_task_numBasisFunctions,          &
                       g09_task_numBasisFunctions ) )

    ! read the overlap matrix 
    do
        read( fno, '(a)', end = 1001 ) read_line
        if ( index( read_line , '*** Overlap ***' ) .ne. 0 ) then
            do 
                ! first read the column numbers. The formatting in the
                ! log file requires a little extra work to read these
                ! values
                col = 0
                read( fno, '(9x,5(i12,2x))', end = 102 ) col(:)
                102 continue
                ! read the info in each row
                do
                    read ( fno, '(i8, 5D14.6)', end=103 ) row,    &
                        (overlap( row, col(i) ), i = 1,           &
                             maxval( col(:) - col(1) ) + 1 )
                    103 continue
                    if ( row == g09_task_numBasisFunctions ) exit
                end do
                if ( maxval( col ) == g09_task_numBasisFunctions ) exit
            end do
            exit
        end if
    end do

    ! close the file
    close( fno )

    ! for convience fill in the upper diagonal too
    do i = 1, g09_task_numBasisFunctions
    do j = i+1, g09_task_numBasisFunctions
        overlap( i, j ) = overlap( j, i )
    end do
    end do

    ! return cleanly
    return

    ! print errors if info not found and abort
    1001 continue
    print'(2a)', ' Overlap matrix not found in file ', logf
    print'(2a)', ' Try running Gausion 09 with IOP(3/33=1)'
    stop
  
end subroutine
