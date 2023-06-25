!**********************************************************************!
!    This module contains variables that will be used throughout       !
!    the g09 helper subroutines
!**********************************************************************!
module g09_commonvar
    implicit none

    ! g09 calculation details
    character(100) :: g09_task_title
    character(10)  :: g09_task_type
    character(30)  :: g09_task_method
    character(30)  :: g09_task_basis

    integer :: g09_task_numAtoms
    integer :: g09_task_numElectrons
    integer :: g09_task_numBasisFunctions
    integer :: g09_task_numBasisFunctionsUsed
    integer :: g09_task_numAOShells
    integer :: g09_task_dof

    real(8) :: g09_task_totalEnergy

    ! structure to hold atom information
    type g09_atom_type
        real(8) x, y, z, mass
        integer atomicNum, basisNum
    end type

    type ( g09_atom_type ), allocatable :: g09_atom(:)

    ! structure to hold atomic orbital basis information
    type g09_AOShell_type
        integer typ, map
    end type

    type( g09_AOShell_type ), allocatable :: g09_AOShell(:)
    
    ! definitions for atomic orbital shells, number of basis
    ! functions for a given shell
    integer, parameter :: g09_AOShellDef(-3:3) = (/7,5,4,1,3,6,10/)

    ! matrices
    real(8), allocatable :: hessian(:,:)
    real(8), allocatable :: moc(:,:)
    real(8), allocatable :: moe(:)
    real(8), allocatable :: overlap(:,:)
    real(8), allocatable :: cic(:,:)
    real(8), allocatable :: fock(:,:)

    ! constants
    ! pi
    real(8), parameter :: pi = 4.d0*datan(1.d0)             
    ! atomic massunits to electron mass units
    real(8), parameter :: AMU_to_electronMass = 1822.88839d0
    ! hartree to wavenumber
    real(8), parameter :: hartree_to_cm = 219474.6313702d0
    ! hartree per boher**2 to mdyne per Angstrom
    real(8), parameter :: au_to_mdyneA  = 15.569141d0
    ! au to debye
    real(8), parameter :: au_to_debye = 2.5415803d0
    ! bohr to angstrom
    real(8), parameter :: bohr_to_angstrom = 0.529177249
                        
end module
