!=============================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
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
!=============================================================================!

module molbuilder_classify
  use crest_parameters,only:wp
  use strucrd,only:coord
  use adjacency
  use canonical_mod
  implicit none
  public

  type,extends(coord) :: coord_classify
    !> new components that are added to the coord type:
    integer,allocatable :: A(:,:)  !> molecular graph/adjacency matrix
    integer,allocatable :: Ah(:,:) !> heavy-atom molecular graph/adjacency

    !> per-atom properties
    real(wp),allocatable :: CN(:)  !> coordination number
    integer,allocatable :: hyb(:) !> hybridization/neighbours count
    integer,allocatable :: nhn(:) !> non-H-neighbours count
    integer,allocatable :: prio(:) !> "invariants"/atom priorities
    logical,allocatable :: inring(:) !> atom part of ring?

  contains
    procedure :: as_coord
    procedure :: from_coord
  end type coord_classify

!==============================================================================!
contains  !> MODULE PROCEDURES START HERE
!==============================================================================!

!> BASIC TYPE PROCEDURES

  function as_coord(this) result(mol)
    class(coord_classify),intent(in) :: this
    type(coord) :: mol

    mol%nat = this%nat
    if (allocated(this%at)) mol%at = this%at
    if (allocated(this%xyz)) mol%xyz = this%xyz

    mol%energy = this%energy
    if (allocated(this%comment)) mol%comment = this%comment
    mol%chrg = this%chrg
    mol%uhf = this%uhf
    mol%nbd = this%nbd
    if (allocated(this%bond)) mol%bond = this%bond
    if (allocated(this%lat)) mol%lat = this%lat
    if (allocated(this%qat)) mol%qat = this%qat
    mol%pdb = this%pdb

  end function as_coord

  subroutine from_coord(this,mol)
    class(coord_classify),intent(inout) :: this
    type(coord),intent(in) :: mol

    this%nat = mol%nat
    if (allocated(mol%at)) this%at = mol%at
    if (allocated(mol%xyz)) this%xyz = mol%xyz

    this%energy = mol%energy
    if (allocated(mol%comment)) this%comment = mol%comment
    this%chrg = mol%chrg
    this%uhf = mol%uhf
    this%nbd = mol%nbd
    if (allocated(mol%bond)) this%bond = mol%bond
    if (allocated(mol%lat)) this%lat = mol%lat
    if (allocated(mol%qat)) this%qat = mol%qat
    this%pdb = mol%pdb
  end subroutine from_coord

!=============================================================================!
!#############################################################################!
!=============================================================================!

!> CLASSIFICATION ROUTINES

  subroutine setup_classify(mol,molc)
    implicit none
    type(coord),intent(in) :: mol
    type(coord_classify),intent(out) :: molc

    real(wp),allocatable :: Bmat(:,:)
    logical,allocatable :: rings(:,:)
    type(canonical_sorter),allocatable :: tmpcan
    integer :: nat
    integer :: ii,jj

    !> Initialize
    call molc%from_coord(mol)
    nat = molc%nat

    !> set up CN, and from that topology
    call mol%cn_to_bond(molc%CN,Bmat,'cov')
    call wbo2adjacency(molc%nat,Bmat,molc%A,0.02_wp)
    deallocate (Bmat)

    !> set up other parameters
    allocate (molc%hyb(nat),source=0)
    allocate (molc%inring(nat),source=.false.)
    call check_rings_min(nat,molc%A,rings)
    do ii = 1,nat
      molc%hyb(ii) = sum(molc%A(:,ii))
      if (any(rings(:,ii))) molc%inring(ii) = .true.
    end do

    allocate (tmpcan)
    call tmpcan%init(mol,invtype='apsp+',heavy=.false.)
    call move_alloc(tmpcan%rank,molc%prio)
    call move_alloc(tmpcan%hadjac,molc%Ah)
    deallocate (tmpcan)

    allocate (molc%nhn(nat),source=0)
    do ii = 1,nat
      molc%nhn(ii) = sum(molc%Ah(:,ii))
    end do

  end subroutine setup_classify

!=============================================================================!
!#############################################################################!
!=============================================================================!
end module molbuilder_classify

