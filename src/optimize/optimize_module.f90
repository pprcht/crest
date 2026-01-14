!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

!> This module wrapps the different optimization algorithms,
!> i.e., this is what can be called for geometry opt.

module optimize_module
  use iso_fortran_env,only:wp => real64
  use crest_parameters
  use crest_calculator
  use strucrd
  use ancopt_module
  use gradientdescent_module
  use rfo_module
  use lbfgs_module
  use optimize_utils
  use thermochem_module
  use hessian_reconstruct
  use newton_raphson_module
  use hr_utils
  !use hessian_tools
  implicit none
  private

  public :: optimize_geometry
  public :: print_opt_data

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine optimize_geometry(mol,molnew,calc,etot,grd,pr,wr,iostatus)
    implicit none
    !> Input
    type(coord)    :: mol
    type(calcdata) :: calc
    logical,intent(in)        :: pr
    logical,intent(in)        :: wr
    !> Output
    type(coord)   :: molnew
    integer,intent(out)       :: iostatus
    real(wp),intent(inout)    :: etot
    real(wp),intent(inout)    :: grd(3,mol%nat)
    real(wp),allocatable :: H_init(:,:),freq(:)
    integer :: nat3
    integer :: io


    iostatus = -1
    !> do NOT overwrite original geometry
    !$omp critical
    molnew%at = mol%at
    molnew%xyz = mol%xyz
    molnew%nat = mol%nat
    !$omp end critical
    nat3 = 3*mol%nat

    !> Check for optimization-individual calculation setup
    if (calc%optnewinit) then
      !$omp critical
      call calc%dealloc_params()
      !$omp end critical
    end if

    !> Check if Hessian Reconstruct is called and initialize the type
    if (calc%do_HR) then
      allocate (calc%chess)
      allocate (H_init(nat3,nat3))
      call calc%chess%alloc(mol%nat,calc%hu_steps,calc%hguess,calc%initialize_hr_type)
    end if

    !> initial singlepoint
    if (calc%do_HR) calc%chess%track_step = .false. !this is not tracked to avoid duplicate
    call engrad(molnew,calc,etot,grd,iostatus)
    if (calc%do_HR) calc%chess%track_step = .true.
    !> optimization
    select case (calc%opt_engine)
    case (0)
      call ancopt(molnew,calc,etot,grd,pr,wr,iostatus)
    case (1)
      !> l-bfgs goes here
      !write(stdout,'(a)') 'L-BFGS currently not implemented'
      !stop
      call lbfgs_optimize(molnew,calc,etot,grd,pr,iostatus)
    case (2)
      !> rfo goes here
      call rfopt(molnew,calc,etot,grd,pr,wr,iostatus)
    case (3)
      !> newton-raphson step goes here, this is a newton step with updated hessians, i.e. quasi Newton
      call newton_raphson(molnew,calc,etot,grd,pr,wr,iostatus)
    case (-1)
      call gradientdescent(molnew,calc,etot,grd,pr,wr,iostatus)
    case default
      write (stdout,'(a)') 'Unknown optimization engine!'
      stop
    end select
    molnew%energy = etot

    if (calc%do_HR .and. iostatus .eq. 0) then !> Hessian construction and post-processing happen here, only do it if geometry relaxation successful
      if (calc%full_HR) then

        write (stdout,*)
        write (stdout,*) "THERMO FROM BFGS" !> This is here for full hessian reconstruct
        write (stdout,*)

        call calc_thermo_from_hess(molnew,calc%chess%H,pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,etot)

      else
        
        call hr_initialize_hessian(calc, molnew%at)
        write(stdout,*)
        write(stdout,*)"THERMO FROM INITIALIZED HESSIAN:"
        write(stdout,*)
        H_init(:,:) = calc%chess%hguess_mat(:,:)
        call calc_thermo_from_hess(molnew,H_init,pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,etot)
        call calc%chess%construct_hessian_bfgs()

        write (stdout,*)
        write (stdout,*) "THERMO FROM RECONSTRUCTED HESSIAN:"
        write (stdout,*)

        call calc_thermo_from_hess(molnew,calc%chess%B,pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,etot)
      end if

      call calc%chess%dealloc()
      deallocate (calc%chess)
    end if

    return
  end subroutine optimize_geometry

!========================================================================================!

  subroutine print_opt_data(calc,ich,natoms,tag)
    implicit none
    type(calcdata) :: calc
    integer,intent(in) :: ich
    integer,intent(in),optional :: natoms
    character(len=*),intent(in),optional :: tag
    integer :: tight,nat
    real(wp) :: ethr,gthr
    character(len=:),allocatable :: ttag
    if (present(tag)) then
      ttag = tag
    else
      ttag = ' '
    end if
    if (present(natoms)) then
      nat = natoms
    else
      nat = 0
    end if

    write (ich,'(a,a)',advance='no') ttag,'Optimization engine: '
    select case (calc%opt_engine)
    case (0)
      write (ich,'(a)') 'ANCOPT'
    case (1)
      write (ich,'(a)') 'L-BFGS'
    case (2)
      write (ich,'(a)') 'RFO'
    case (-1)
      write (ich,'(a)') 'Gradient Descent'
    case default
      write (ich,'(a)') 'Unknown'
    end select
    if (calc%opt_engine >= 0) then
      write (ich,'(a,a)',advance='no') ttag,'Hessian update type: '
      select case (calc%iupdat)
      case (0)
        write (ich,'(a)') 'BFGS'
      case (1)
        write (ich,'(a)') 'Powell'
      case (2)
        write (ich,'(a)') 'SR1'
      case (3)
        write (ich,'(a)') 'Bofill'
      case (4)
        write (ich,'(a)') 'Farkas-Schlegel'
      end select
    end if

    tight = calc%optlev
    call get_optthr(nat,tight,calc,ethr,gthr)
    write (ich,'(a,a,e10.3,a,e10.3,a)') ttag,'E/G convergence criteria: ',&
    & ethr,' Eh,',gthr,' Eh/a0'

    write (ich,'(a,a,i0)') ttag,'maximum optimization steps: ',calc%maxcycle

  end subroutine print_opt_data

!========================================================================================!
!========================================================================================!
end module optimize_module
