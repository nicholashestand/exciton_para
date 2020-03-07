!**********************************************************************!
!   Get information about a Gaussian 09 calculation from the           !
!   formatted checkpoint file.                                         !
!**********************************************************************!
subroutine g09_hessian( fch )
    use g09_commonvar
    implicit none

    character(*), intent(in) :: fch
    integer, parameter :: fno = 509 
    logical exists
    character(100) read_line
    integer i, j

    ! check that the formatted checkpoint file exists
    inquire( file = fch, exist = exists )
    if ( .not. exists ) then
        print'(3a)', ' Formatted checkpoint file ', fch, &
                     ' not found.'
        stop
    end if

    ! open the formatted checkpoint file and read the hessian matrix
    open( unit = fno, file = fch, status = 'old', &
          action = 'read' )

    ! allocate space for the hessian
    allocate( hessian( g09_task_dof, g09_task_dof) )

    ! read the hessian matrix in upper triangular form
    do
        read( fno, '(a)', end = 1001 ) read_line
        if ( index( read_line , 'Cartesian Force Constants' ) .ne. 0 ) then
            read( fno, '(5E16.8)' ) ((hessian(j,i), j = 1, i  ), i = 1, g09_task_dof )
            exit
        end if
    end do

    ! close the file
    close( fno )

    ! for convience fill in the lower diagonal too
    do i = 1, g09_task_dof
    do j = 1, i-1
        hessian(i,j) = hessian(j,i)
    end do
    end do

    ! return cleanly
    return

    ! print errors if info not found and abort
    1001 continue
    print'(2a)', 'Cartesian force constants not found in file ', fch
    stop
  
end subroutine
