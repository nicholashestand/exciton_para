!**********************************************************************!
!   Calculate the charge transfer integrals from the output of a       !
!   Gaussian calculation.                                              !
!**********************************************************************!
program g09_ctint
    use g09_commonvar
    implicit none

    character(200) :: fch_m1, log_m1, fch_m2, log_m2,   &
                      fch_d, log_d, fout
    real(8), allocatable :: s1(:,:), s2(:,:), sd(:,:),          &
                            m1moc(:,:), m2moc(:,:), mdmoc(:,:), &
                            m12moc(:,:), mdmoe(:,:), SMO(:,:),  & 
                            mdmocinv(:,:), HAO(:,:), HMO(:,:)
    real(8) Hsub(4,4), Ssub(4,4)
    integer homo1, homo2, lumo1, lumo2, i, baseN1, baseN2, baseNd, j
    integer, parameter :: fno = 103


    ! initialize the program, read options and input file
    call ctint_init(fch_m1, log_m1, fch_m2, log_m2, &
                    fch_d, log_d, fout)
    
    ! read in the data from the Guassian files
    print*, '>> Getting info about molecule 1.'
    call g09_calcInfo(fch_m1)
    ! alpha and beta spin electrons are the same in the G09 calcs and
    ! only one set of orbitals is calculated for the setup given with
    ! the -m option
    baseN1 = g09_task_numBasisFunctions
    homo1 = ceiling(g09_task_numElectrons/2.d0)
    lumo1 = homo1 + 1
    call g09_mocoeff(fch_m1)
    call g09_overlap(log_m1)
    ! allocate space for the matrices
    allocate( s1( baseN1, baseN1 ), m1moc( baseN1, baseN1 ))
    s1 = overlap
    m1moc = moc

    print*, '>> Getting info about molecule 2.'
    call g09_calcInfo(fch_m2)
    baseN2 = g09_task_numBasisFunctions
    homo2 = ceiling(g09_task_numElectrons/2.d0) + baseN1
    lumo2 = homo2 + 1
    call g09_mocoeff(fch_m2)
    call g09_overlap(log_m2)
    ! allocate space for the matrices
    allocate( s2( baseN2, baseN2 ), m2moc( baseN2, baseN2 ))
    s2 = overlap
    m2moc = moc

    print*, '>> Getting info about the dimer.'
    call g09_calcInfo(fch_d)
    baseNd = g09_task_numBasisFunctions
    call g09_mocoeff(fch_d)
    call g09_overlap(log_d)
    ! allocate space for the matrices
    allocate( sd( baseNd, baseNd ) , mdmoc( baseNd, baseNd ),&
           m12moc( baseNd, baseNd ), mdmoe( baseNd, baseNd ),&
           HAO( baseNd, baseNd ), HMO( baseNd, baseNd ),     &
           SMO( baseNd, baseNd ), mdmocinv( baseNd, baseNd ) )
    sd = overlap
    mdmoc = moc
    mdmoe = 0.d0
    do i = 1, baseNd
        mdmoe(i,i) = moe(i)*Hartree_to_cm
    end do

    ! combine the matrices for the monomer mos into one matrix
    m12moc = 0.d0
    m12moc( 1:baseN1, 1:baseN1 ) = m1moc(:,:)
    m12moc( baseN1+1:baseN1+baseN2, baseN1+1:baseN1+baseN2 ) = m2moc(:,:)

    ! calculate the dimer hamiltonian in the atomic orbital basis
    ! by back transforming the eigenvalues
    call find_inverse( mdmoc, baseNd, mdmocinv )
    HAO = matmul( matmul( matmul( sd, mdmoc ) , mdmoe ) , mdmocinv )

    ! mow project the Hamiltonian onto the molecular orbital basis set
    HMO = matmul( matmul( transpose(m12moc), HAO ), m12moc )
    ! and the overlap matrix
    SMO = matmul( matmul( transpose(m12moc), sd ), m12moc )

    ! get submatrices with just the important MOs for us
    call get_submatrix( HMO, SMO, Hsub, Ssub, homo1, homo2, lumo1, lumo2, &
                        g09_task_numBasisFunctions )

    ! write out the orthogonalized and nonorthogonalized Hamiltonian
    open( unit = fno, file = trim(fout), action='write')
    write( fno, * ) 'Nonorthogonal Sub Hamiltonian (cm-1)'
    write( fno, * ) ',|HOMO1>,|HOMO2>,|LUMO1>,|LUMO2>'
    write( fno, '(a,4(",",f14.4))' ) '<HOMO1|',(Hsub(1,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<HOMO2|',(Hsub(2,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<LOMO1|',(Hsub(3,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<LOMO2|',(Hsub(4,i),i=1,4)

    call orthogonalize( Hsub, Ssub, 4 )

    write( fno, * ) 'Orthogonalized Sub Hamiltonian (cm-1)'
    write( fno, * ) ',|HOMO1>,|HOMO2>,|LUMO1>,|LUMO2>'
    write( fno, '(a,4(",",f14.4))' ) '<HOMO1|',(Hsub(1,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<HOMO2|',(Hsub(2,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<LOMO1|',(Hsub(3,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<LOMO2|',(Hsub(4,i),i=1,4)

    write( fno, * ) 'Overlap Sub Matrix'
    write( fno, * ) ',|HOMO1>,|HOMO2>,|LUMO1>,|LUMO2>'
    write( fno, '(a,4(",",f14.4))' ) '<HOMO1|',(Ssub(1,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<HOMO2|',(Ssub(2,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<LOMO1|',(Ssub(3,i),i=1,4)
    write( fno, '(a,4(",",f14.4))' ) '<LOMO2|',(Ssub(4,i),i=1,4)

    close( fno )

    print*, 
    print*, '>> Wrote output to file ', trim(adjustl(fout))
    print*, 

end program


!**********************************************************************!
!   orthogonalize a matrix H given its overlap matrix S                !
!**********************************************************************!
subroutine orthogonalize( H, S, N )
    implicit none

    integer, intent(in) :: N
    real(8), intent(in) :: S(N,N)
    real(8), intent(inout) :: H(N,N)
    real(8) eval(N), evalMat(N,N), SWork(N,N)
    integer i

    ! take square root of overlap matrix
    Swork = S
    call diagonalize( SWork, N, eval )
    evalMat = 0.d0
    do i = 1, N
        evalMat(i,i) = 1.d0/dsqrt(eval(i))
    end do
    ! convert back to the original basis
    SWork = matmul( SWork, matmul( evalMat, transpose( SWork ) ) )

    ! Now orthogonalize the matrix H using the lowdin method
    H = matmul( SWork, matmul(H, SWork ) )

end subroutine


!**********************************************************************!
!   get the submatrices of the orbitals of interest                    !
!**********************************************************************!
subroutine get_submatrix( HMO, SMO, Hsub, Ssub, homo1, homo2, &
                          lumo1, lumo2, d )
    implicit none

    integer, intent(in) :: d, homo1, homo2, lumo1, lumo2
    real(8), intent(in) :: HMO(d,d), SMO(d,d)
    real(8), intent(out) :: Hsub(4,4), Ssub(4,4)
    integer i, j

    ! get the subspace of the Hamiltonian corresponding
    ! to the HOMOs and LUMOs. 
    Hsub( 1, 1 ) = HMO( homo1, homo1 )
    Hsub( 1, 2 ) = HMO( homo1, homo2 )
    Hsub( 1, 3 ) = HMO( homo1, lumo1 )
    Hsub( 1, 4 ) = HMO( homo1, lumo2 )
    Hsub( 2, 2 ) = HMO( homo2, homo2 )
    Hsub( 2, 3 ) = HMO( homo2, lumo1 )
    Hsub( 2, 4 ) = HMO( homo2, lumo2 )
    Hsub( 3, 3 ) = HMO( lumo1, lumo1 )
    Hsub( 3, 4 ) = HMO( lumo1, lumo2 )
    Hsub( 4, 4 ) = HMO( lumo2, lumo2 )

    Ssub( 1, 1 ) = SMO( homo1, homo1 )
    Ssub( 1, 2 ) = SMO( homo1, homo2 )
    Ssub( 1, 3 ) = SMO( homo1, lumo1 )
    Ssub( 1, 4 ) = SMO( homo1, lumo2 )
    Ssub( 2, 2 ) = SMO( homo2, homo2 )
    Ssub( 2, 3 ) = SMO( homo2, lumo1 )
    Ssub( 2, 4 ) = SMO( homo2, lumo2 )
    Ssub( 3, 3 ) = SMO( lumo1, lumo1 )
    Ssub( 3, 4 ) = SMO( lumo1, lumo2 )
    Ssub( 4, 4 ) = SMO( lumo2, lumo2 )

    do i = 1, 4
    do j = i+1, 4
        Hsub( j, i ) = Hsub(i,j)
        Ssub( j, i ) = Ssub(i,j)
    end do
    end do

end subroutine



!**********************************************************************!
!   Invert a matrix using lapack                                       !
!**********************************************************************!
subroutine find_inverse( A, N, W )
    implicit none

    integer, intent(in) :: N
    real(8), intent(in) :: A(N,N)
    real(8), intent(out):: W(N,N)

    real*8 work( N, N )
    integer ipiv( N ), info

    W = A

    call dgetrf( N, N, W, N, ipiv, info )
    call dgetri( N, W, N, ipiv, work, N, info )

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
!   Initialize the program. Read command line options and input file   ! 
!**********************************************************************!
subroutine ctint_init(fch_m1, log_m1, fch_m2, log_m2, fch_d, log_d, &
                      fout )
    implicit none

    character(200), intent(out) :: fch_m1, log_m1, fch_m2, log_m2, &
                                   fch_d, log_d, fout
    integer nargs, narg, ios, line, pos
    integer, parameter :: fno = 67, fno2 = 68
    character(2) nproc
    character(100) arg, fin, label, fxyzm1, fxyzm2, task, method
    character(100) buff
    logical exists, makeinput

    makeinput = .false.
    fin = ''
    fout = ''
    fxyzm1 = ''
    fxyzm2 = ''
    nproc = '1'

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
                    case('task')
                        read( buff, *, iostat = ios ) task
                        print*, 'task : ', trim(adjustl(task))
                    case('fch_m1')
                        read( buff, *, iostat = ios ) fch_m1
                        print*, 'fch_m1: ', trim(adjustl(fch_m1))
                    case('log_m1')
                        read( buff, *, iostat = ios ) log_m1
                        print*, 'log_m1: ', trim(adjustl(log_m1))
                    case('fch_m2')
                        read( buff, *, iostat = ios ) fch_m2
                        print*, 'fch_m2: ', trim(adjustl(fch_m2))
                    case('log_m2')
                        read( buff, *, iostat = ios ) log_m2
                        print*, 'log_m2: ', trim(adjustl(log_m2))
                    case('fch_d')
                        read( buff, *, iostat = ios ) fch_d
                        print*, 'fch_d: ', trim(adjustl(fch_d))
                    case('log_d')
                        read( buff, *, iostat = ios ) log_d
                        print*, 'log_d: ', trim(adjustl(log_d))
                    case('xyz_file_m1')
                        read( buff, *, iostat = ios ) fxyzm1
                        print*, 'fxyzm1: ', trim(adjustl(fxyzm1))
                    case('xyz_file_m2')
                        read( buff, *, iostat = ios ) fxyzm2
                        print*, 'fxyzm1: ', trim(adjustl(fxyzm2))
                    case('nproc')
                        read( buff, *, iostat = ios ) nproc
                        print*, 'nproc : ', trim(adjustl(nproc))
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
        if ( method == '' ) method = 'HF/6-31G'
        
        ! make sure the xyz files exists
        inquire( file = trim(fxyzm1), exist = exists )
        if ( .not. exists ) then
            print*, 'The file ', trim(adjustl(fxyzm1)), &
                    ' does not exist. Aborting...'
            stop
        end if

        inquire( file = trim(fxyzm2), exist = exists )
        if ( .not. exists ) then
            print*, 'The file ', trim(adjustl(fxyzm2)), &
                    ' does not exist. Aborting...'
            stop
        end if
 
        ! the calculation for molecule 1
        open( unit = fno, file = trim(adjustl(task))//'_m1.gjf', &
              action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))//'_m1.chk'
        write( fno, * ) '%NProcShared='//trim(adjustl(nproc))
        write( fno, * ) '%mem=1GB'
        write( fno, * ) '# SP '//trim(method)//' IOP(3/33=1) IOP(6/7=3) '//&
                        'Symmetry=None MaxDisk=2GB'
        write( fno, * )
        write( fno, * ) trim(adjustl(task))//'_m1'
        write( fno, * )
        write( fno, * ) '0 1'
        ! write the coordinates from the xyz file
        open( unit = fno2, file = trim(fxyzm1), action = 'read' )
        read( fno2, * )
        read( fno2, * )
        do
            read( fno2, '(a)', end=105) buff
            write( fno, * ) buff
        end do
105     continue
        ! write a blank line so g09 doesnt crash
        write( fno, * ) 
        close( fno2 )
        close( fno )

        ! the calculation for molecule 2
        open( unit = fno, file = trim(adjustl(task))//'_m2.gjf', &
              action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))//'_m2.chk'
        write( fno, * ) '%NProcShared='//trim(adjustl(nproc))
        write( fno, * ) '%mem=1GB'
        write( fno, * ) '# SP '//trim(method)//' IOP(3/33=1) IOP(6/7=3) '//&
                        'Symmetry=None MaxDisk=2GB'
        write( fno, * )
        write( fno, * ) trim(adjustl(task))//'_m2'
        write( fno, * )
        write( fno, * ) '0 1'
        ! write the coordinates from the xyz file
        open( unit = fno2, file = trim(fxyzm2), action = 'read' )
        read( fno2, * )
        read( fno2, * )
        do
            read( fno2, '(a)', end=106) buff
            write( fno, * ) buff
        end do
106     continue
        ! write a blank line so g09 doesnt crash
        write( fno, * ) 
        close( fno2 )
        close( fno )

        ! the calculation for the dimer
        open( unit = fno, file = trim(adjustl(task))//'_d.gjf', &
              action = 'write' )
        write( fno, * ) '%Chk='//trim(adjustl(task))//'_d.chk'
        write( fno, * ) '%NProcShared='//trim(adjustl(nproc))
        write( fno, * ) '%mem=1GB'
        write( fno, * ) '# SP '//trim(method)//' IOP(3/33=1) IOP(6/7=3) '//&
                        'Symmetry=None MaxDisk=2GB'
        write( fno, * )
        write( fno, * ) trim(adjustl(task))//'_d'
        write( fno, * )
        write( fno, * ) '0 1'
        ! write the coordinates from the xyz file
        open( unit = fno2, file = trim(fxyzm1), action = 'read' )
        read( fno2, * )
        read( fno2, * )
        do
            read( fno2, '(a)', end=107) buff
            write( fno, * ) buff
        end do
107     continue
        close( fno2 )
        open( unit = fno2, file = trim(fxyzm2), action = 'read' )
        read( fno2, * )
        read( fno2, * )
        do
            read( fno2, '(a)', end=108) buff
            write( fno, * ) buff
        end do
108     continue
        ! write a blank line so g09 doesnt crash
        write( fno, * ) 
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
    print*, ' g09_ctint version 1.0'
    print*, 
    print*, ' This  program  calculates  the charge  transfer  integrals '
    print*, ' between two molecules by transforming the Hamiltonian into '
    print*, ' the molecular orbital basis.  The calculation  uses output '
    print*, ' from a Gaussian 09 calculation.'
    print*,
    print*, ' Usage:'
    print*,
    print*, ' -- help, -h: print this message'
    print*, ' -i file    : specifies the input file containing the names '
    print*, '              of the formatted checkpoint files. If the  -m '
    print*, '              option is also given, it reads the input file '
    print*, '              for information about the G09 calculation  to '
    print*, '              add to the Gaussian 09 input files.'
    print*, ' -o file    : specifies the output file.'
    print*, ' -m         : makes three gaussian input files from two xyz '
    print*, '              files of atomic coordinates. The three  files '
    print*, '              are  for  two  monomer  calculations  and one '
    print*, '              dimer calculation. The  options  are  set  so '
    print*, '              that the resulting fch and output  files  can '
    print*, '              be used directly with this program. Note that '
    print*, '              the formatted checkpoint  file  needs  to  be '
    print*, '              made by the Gaussian 09 formchk utility.'
    print*, 
    print*, ' To  use  this  program  to  calculate  the charge transfer '
    print*, ' integrals, you must first run single point calculations in '
    print*, ' Gaussian 09 of the two molecules of interest and  a single '
    print*, ' point calculation of the dimer system. The  -m  option  of '
    print*, ' this program can  be  used  to  generate  the  appropriate '
    print*, ' Gaussian 09 input files if xyz  coordinate files  for  the '
    print*, ' molecules of interest are provided.'
    print*, 
    print*, ' The  resulting  checkpoint  files  must  be  converted  to '
    print*, ' formatted checkpoint files using the Gaussian 09 formcheck '
    print*, ' utility before they can be read by this program.'
    print*,

end subroutine


