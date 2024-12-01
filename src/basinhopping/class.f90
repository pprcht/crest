!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2024 Philipp Pracht, David Wales
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

module bh_class_module
  use crest_parameters
  use strucrd,only:coord
  implicit none

!=========================================================================================!

  public :: bh_class
  type :: bh_class
!************************************************************************
!* data object that contains the data for a *SINGLE* basin-hopping chain
!************************************************************************

!>--- counters
    integer :: iteration = 0   !> current iteration
    integer :: saved = 0       !> number of saved quenches

!>--- paramters
    integer  :: maxsteps = 100     !> maximum steps to take
    real(wp) :: temp = 300.0_wp    !> MC acceptance temperature
    real(wp) :: scalefac = 1.0_wp  !> temperature increase factor
    integer :: steptype = 0        !> step type selection
    real(wp) :: stepsize(3) = &    !> step sizes e.g. for lengths, angles, dihedrals
    &    (/0.2_wp,0.2_wp,0.2_wp/)
    integer :: maxsave = 100       !> maximum number of quenches saved

!>--- results/properties
    real(wp) :: emin = 0.0_wp  !> current ref energy of markov chain
    integer  :: whichmin = 0   !> mapping to which structure emin refers
    real(wp) :: emax = 0.0_wp  !> highest energy structure among saved quenches
    integer  :: whichmax = 0   !> mapping of highest energy structure
    type(coord),allocatable :: structures(:)  !> list of structures from succesfull quenches

!>--- Type procedures
  contains
    procedure :: init => bh_class_allocate
    procedure :: deallocate => bh_class_deallocate
    procedure :: add => bh_class_add
  end type bh_class

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine bh_class_allocate(self,temp,maxsteps,maxsave)
    implicit none
    class(bh_class) :: self
    real(wp),intent(in),optional :: temp
    integer,intent(in),optional  :: maxsteps
    integer,intent(in),optional  :: maxsave

    call self%deallocate()
    if (present(temp)) then
      self%temp = temp
    end if
    if (present(maxsteps)) then
      self%maxsteps = maxsteps
    end if
    if (present(maxsave)) then
      self%maxsave = maxsave
    end if

    self%iteration = 0
    self%saved = 0
    allocate (self%structures(self%maxsave))
  end subroutine bh_class_allocate

!=========================================================================================!

  subroutine bh_class_deallocate(self)
    implicit none
    class(bh_class) :: self
    if (allocated(self%structures)) deallocate (self%structures)
  end subroutine bh_class_deallocate

!=========================================================================================!

  subroutine bh_class_add(self,mol)
    implicit none
    class(bh_class) :: self
    type(coord) :: mol
    integer :: i,j
    real(wp) :: mintmp,maxtmp
    if(self%saved < self%maxsave)then
      self%saved = self%saved + 1
      self%structures( self%saved ) = mol
    else
      i = self%%whichmax
      self%structures( i ) = mol
    endif

    mintmp = huge(mintmp)
    maxtmp = -huge(maxtmp)
    do i = 1,self%saved
      if(structures(i)%energy < mintmp)then
        mintmp = structures(i)%energy
        self%whichmin = i
      endif
      if(structures(i)%energy > maxtmp)then
        maxtmp = structures(i)%energy
        self%whichmax = i
      endif
    enddo
    self%emin = mintmp
    self%emax = maxtmp
  end subroutine bh_class_add

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!
end module bh_class_module
