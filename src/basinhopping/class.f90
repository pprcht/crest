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
  use canonical_mod
  use irmsd_module
  implicit none

!=========================================================================================!

  public :: bh_class
  type :: bh_class
!************************************************************************
!* data object that contains the data for a *SINGLE* basin-hopping chain
!************************************************************************
    integer :: id = 0 !> Run/Thread ID

!>--- counters
    integer :: iteration = 0   !> current iteration
    integer :: saved = 0       !> number of saved quenches

!>--- paramters
    integer  :: maxsteps = 100     !> maximum steps to take
    real(wp) :: temp = 300.0_wp    !> MC acceptance temperature
    real(wp) :: scalefac = 1.0_wp  !> temperature increase factor
    real(wp) :: rthr = 0.125_wp    !> RMSD threshold (\AA)
    real(wp) :: ethr = 0.05_wp     !> minima/conformer energy distinction (kcal/mol)
    integer :: steptype = 0        !> step type selection
    real(wp) :: stepsize(3) = &    !> step sizes e.g. for lengths, angles, dihedrals
    &    (/0.2_wp,0.2_wp,0.2_wp/)
    integer :: maxsave = 100       !> maximum number of quenches saved
    real(wp),allocatable :: etarget  !> target energy to be hit (useful in benchmarks)

!>--- results/properties
    real(wp) :: emin = 0.0_wp  !> current ref energy of markov chain
    integer  :: whichmin = 0   !> mapping to which structure emin refers
    real(wp) :: emax = 0.0_wp  !> highest energy structure among saved quenches
    integer  :: whichmax = 0   !> mapping of highest energy structure
    type(coord),allocatable :: structures(:)  !> list of structures from succesfull quenches

!>--- temporary storage
    integer,allocatable  :: amat(:,:)        !> adjacency matrix
    real(wp),allocatable :: zmat(:,:)        !> internal coordinates (to cache the memory)
    type(rmsd_cache),allocatable :: rcache   !> similarity check cache (iRMSD)
    logical :: stereocheck = .false.         !> check for false-rotamers?
    type(canonical_sorter),allocatable :: sorters(:) !> canonical atom ID storage
    logical :: topocheck = .true.            !> check for correct connectivity
    type(canonical_sorter),allocatable :: refsort  !> use same reference connectivity for all

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
    self%maxsave = min(self%maxsave,self%maxsteps)

    self%iteration = 0
    self%saved = 0
    allocate (self%structures(self%maxsave))
    allocate (self%sorters(self%maxsave))
  end subroutine bh_class_allocate

!=========================================================================================!

  subroutine bh_class_deallocate(self)
    implicit none
    class(bh_class) :: self
    if (allocated(self%structures)) deallocate (self%structures)
    if (allocated(self%amat)) deallocate (self%amat)
    if (allocated(self%zmat)) deallocate (self%zmat)
    if (allocated(self%sorters)) deallocate (self%sorters)
    if (allocated(self%rcache)) deallocate(self%rcache)
    if (allocated(self%refsort)) deallocate(self%refsort)
  end subroutine bh_class_deallocate

!=========================================================================================!

  subroutine bh_class_add(self,mol)
    implicit none
    class(bh_class) :: self
    type(coord) :: mol
    integer :: i,j
    real(wp) :: mintmp,maxtmp
    if (self%saved < self%maxsave) then
      self%saved = self%saved+1
      i = self%saved
      self%structures(i) = mol
      !$omp critical
      call self%sorters(i)%init(mol,invtype='apsp+',heavy=.false.)
      if (i == 1) then
        self%stereocheck = .not. (self%sorters(i)%hasstereo(mol))
      end if
      !$omp end critical
    else
      i = self%whichmax
      self%structures(i) = mol
      !$omp critical
      call self%sorters(i)%deallocate()
      call self%sorters(i)%init(mol,invtype='apsp+',heavy=.false.)
      call self%sorters(i)%shrink()
      !$omp end critical
    end if

    mintmp = huge(mintmp)
    maxtmp = -huge(maxtmp)
    do i = 1,self%saved
      if (self%structures(i)%energy < mintmp) then
        mintmp = self%structures(i)%energy
        self%whichmin = i
      end if
      if (self%structures(i)%energy > maxtmp) then
        maxtmp = self%structures(i)%energy
        self%whichmax = i
      end if
    end do
    self%emin = mintmp
    self%emax = maxtmp
  end subroutine bh_class_add

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!
end module bh_class_module
