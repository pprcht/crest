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
!================================================================================!

!====================================================!
! module lammps_interface
! An interface to the LAMMPS calculator
!====================================================!

module lammps_interface
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
#ifdef WITH_LAMMPS
   USE LIBLAMMPS
#endif
  implicit none
  private

#ifndef WITH_LAMMPS
  !> these are placeholders if no lammps module is used!
  type :: lammps
    integer :: id = 0
  end type lammps
#endif

  public :: lammps
  public :: lmp_setup
  public :: lmp_sp
  public :: lmp_getwbos
  public :: lmp_print


  logical,parameter,private :: DEBUG = .true. 

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine lmp_setup(mol,lmp,infile)
!********************************************
!* subroutine lmp_setup
!* Initialize the lammps calculator
!* Since this involves some I/O, it must be
!* called in a threadsafe region of the code
!********************************************
    implicit none
    type(coord),intent(in)  :: mol
    type(lammps),intent(inout) :: lmp
    character(len=*),intent(in) :: infile
#ifdef WITH_LAMMPS
!>--- create a new lmp object
    lmp = lammps()
    if(DEBUG)then 
      PRINT*, 'LAMMPS Version: ', lmp%version()
    endif
!>--- read the input file
    CALL lmp%file(infile)

     

   
#else /* WITH_LAMMPS */
     call lmp_avail()
#endif
  end subroutine lmp_setup

!========================================================================================!

  subroutine lmp_sp(mol,lmp,energy,gradient,iostatus)
    implicit none
    !> INPUT
    type(coord),intent(in)  :: mol
    type(lammps),intent(inout) :: lmp
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,mol%nat)
    integer,intent(out) :: iostatus
    !> LOCAL
    logical :: fail 
    energy = 0.0_wp
    gradient = 0.0_wp
    iostatus = 0
    fail = .false.
#ifdef WITH_LAMMPS
!>--- update the coordinates for the lmp object
    call lmp%scatter_atoms('x', reshape(mol%xyz, [3*mol%nat]))    

!>--- run a single evaluation
    call lmp%command('run 0')

!>--- read the computed potential energy
    energy = lmp%get_thermo('pe')

!>--- read the calculated forces from the lmp type
    call lmp%gather_atoms('f',3,gradient)
    gradient = -gradient

    if(fail)then
      iostatus = -1
    endif
#else
    call lmp_avail()
#endif
  end subroutine lmp_sp

!========================================================================================!

   subroutine lmp_print(iunit,lmp)
    implicit none
    integer,intent(in) :: iunit
    type(lammps),intent(in) :: lmp
#ifdef WITH_LAMMPS

#endif
    return
   end subroutine lmp_print


!========================================================================================!
!> obtain bond orders from lammps
  subroutine lmp_getwbos(lmp,nat,wbo)
    implicit none
    type(lammps),intent(in) :: lmp
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)
    integer ndim
    wbo = 0.0_wp
#ifdef WITH_LAMMPS

#endif   
  end subroutine lmp_getwbos

!========================================================================================!

  subroutine lmp_avail()
     implicit none
#ifndef WITH_LAMMPS
    write (stdout,*) 'Error: Compiled without LAMMPS support!'
    write (stdout,*) 'Use -DWITH_LAMMPS=true in the setup to enable this function'
    error stop
#endif
  end subroutine lmp_avail

!========================================================================================!
!========================================================================================!
end module lammps_interface

