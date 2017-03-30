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

end module
