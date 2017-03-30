!**********************************************************************!
!   Read in the CI matrix from a Gaussian output file                  !
!   The excited state of interest must be passed to this subroutine    !
!**********************************************************************!
subroutine g09_cicoeff( logf, exstate )
    use g09_commonvar
    implicit none

    character(*), intent(in) :: logf
    integer, intent(in) :: exstate
    integer, parameter :: fno = 509 
    logical exists
    character(100) read_line, matchline
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
    allocate( cic( g09_task_numBasisFunctions,          &
                   g09_task_numBasisFunctions ) )

    ! read the ci matrix 
    cic = 0.d0
    ! write the line to search for in the log file
    write( matchline, '(a13,i4)' ) 'Excited State', exstate
    do
        read( fno, '(a)', end = 1001 ), read_line
        if ( index( read_line , trim(matchline) ) .ne. 0 ) then
            do 
                read( fno, '(a)' ), read_line
                ! store ci coefficients, or exit if we have read past them
                if ( index( read_line, '->' ) .ne. 0 ) then
                    read( read_line, '(i9,2x,i9,f10.5)' ) i, j, cic( i, j )
                else if ( index( read_line, '<-') .ne. 0 ) then
                    read( read_line, '(i9,2x,i9,f10.5)' ) i, j, cic( j, i )
                else
                    exit
                end if
            end do
            exit
        end if
    end do

    ! normalize the ci coefficients to 1, the are normalized to 0.5 by
    ! default in unrestricted gaussian calculations
    cic = cic * dsqrt(2.d0)

    ! close the file
    close( fno )

    ! return cleanly
    return

    ! print errors if info not found and abort
    1001 continue
    print'(2a)', ' CI Coefficients not found in file ', logf
    print'(2a)', ' Try running a CIS calculation in Gausion 09'
    stop
  
end subroutine
