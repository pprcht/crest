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

module molbuilder_classify_type
  use crest_parameters,only:wp,stdout
  use strucrd,only:coord,i2e,sumform
  use adjacency
  use canonical_mod
  implicit none
  private

  type :: functional_group
    character(len=:),allocatable :: name
    integer :: natms = 0
    integer,allocatable :: ids(:)
    integer :: attached_to = 0
    logical :: seeded = .false.
  contains
    procedure :: clear => clear_func_group
    procedure :: copy => copy_func_group
  end type functional_group

  type,extends(coord) :: coord_classify
    !> new components that are added to the coord type:
    integer,allocatable :: A(:,:)  !> molecular graph/adjacency matrix
    integer,allocatable :: Ah(:,:) !> heavy-atom molecular graph/adjacency

    !> per-atom properties/information
    real(wp),allocatable :: CN(:)               !> coordination number
    integer,allocatable :: hyb(:)               !> hybridization/neighbours count
    integer,allocatable :: nhn(:)               !> non-H-neighbours count
    integer,allocatable :: prio(:)              !> "invariants"/atom priorities
    logical,allocatable :: inring(:)            !> atom part of ring?
    logical,allocatable :: term(:)              !> terminal atom (H,F,Cl,...,=O,etc.)
    character(len=10),allocatable :: atinfo(:)  !> atom info

    !> functional groups
    integer :: nfuncs=0
    type(functional_group),allocatable :: funcgroups(:)

    !> utility storage
    logical,allocatable :: lwork(:)
    integer,allocatable :: iwork(:)

  contains
    procedure :: as_coord
    procedure :: from_coord
    generic,public :: add => coord_classify_add_fg
    procedure,private :: coord_classify_add_fg
    procedure :: print_funcgroups => coord_classify_print_functional
  end type coord_classify

  public :: coord_classify   !> the extended coord type
  public :: functional_group !> subtype of coord_classify
  public :: setup_classify   !> setup a coord_classify from coord
  public :: atinfo_classify  !> add atinfo string to a coord_classify

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

  subroutine clear_func_group(self)
    implicit none
    class(functional_group) :: self
    if (allocated(self%name)) deallocate (self%name)
    if (allocated(self%ids)) deallocate (self%ids)
    self%attached_to = 0
    self%natms = 0
    self%seeded = .false.
  end subroutine clear_func_group

  subroutine copy_func_group(self,fg)
    implicit none
    class(functional_group) :: self
    type(functional_group) :: fg
    if (allocated(fg%name)) self%name = fg%name
    if (allocated(fg%ids)) self%ids = fg%ids
    self%attached_to = fg%attached_to
    self%natms = fg%natms
    self%seeded = fg%seeded
  end subroutine copy_func_group

  subroutine coord_classify_add_fg(self,fg)
    implicit none
    class(coord_classify) :: self
    type(functional_group) :: fg
    type(functional_group),allocatable :: fg_list(:)
    integer :: ii,jj
    if (.not.allocated(self%funcgroups)) then
      allocate (self%funcgroups(1))
      call self%funcgroups(1)%copy(fg)
    else
      ii = size(self%funcgroups,1)
      allocate (fg_list(ii+1))
      do jj = 1,ii
        call fg_list(jj)%copy(self%funcgroups(jj))
      end do
      call fg_list(ii+1)%copy(fg)
      call move_alloc(fg_list,self%funcgroups)
    end if
    self%nfuncs = size(self%funcgroups,1)
  end subroutine coord_classify_add_fg

!=============================================================================!
!#############################################################################!
!=============================================================================!

!> CLASSIFICATION ROUTINES

  subroutine setup_classify(mol,molc)
    !***************************************************
    !* set up the derived coord_classify object "molc"
    !* from a standard coord object "mol".
    !* in particular, adjacency graphs, CN, etc.
    !***************************************************
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
    allocate (molc%term(nat),source=.false.)
    call check_rings_min(nat,molc%A,rings)
    do ii = 1,nat
      molc%hyb(ii) = sum(molc%A(:,ii))
      if (any(rings(:,ii))) molc%inring(ii) = .true.
      if (molc%hyb(ii) .eq. 1) molc%term(ii) = .true.
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

  subroutine atinfo_classify(molc)
    !*****************************************
    !* Update a coord_classify object "molc"
    !* and fill in its atinfo strings based
    !* on some basic chemoinformatics.
    !****************************************
    implicit none
    type(coord_classify),intent(inout) :: molc
    integer :: nat
    integer :: ii,jj

    if (molc%nat <= 0) then
      write (stdout,*) 'molc not allocated in atinfo_classify()'
      return
    end if
    if (allocated(molc%atinfo)) deallocate (molc%atinfo)
    nat = molc%nat
    allocate (molc%atinfo(nat),source=repeat(' ',10))

    do ii = 1,molc%nat
      associate (str => molc%atinfo(ii),ati => molc%at(ii))
        str = trim(i2e(ati))
        select case (ati)

        case (6) !> carbon
          if (molc%hyb(ii) == 3) then !> sp2
            str = trim(i2e(ati,'lc'))
          else if (molc%hyb(ii) == 4) then !> sp3
            if (molc%nhn(ii) == 1) then
              str = trim(str)//'H3'
            else if (molc%nhn(ii) == 2) then
              str = trim(str)//'H2'
            else if (molc%nhn(ii) == 0) then
              str = 'methane'
            end if
          end if

        case (7) !> nitrogen
          if (molc%hyb(ii) == 3) then
            if (molc%nhn(ii) == 1) then
              str = trim(str)//'H2'
            else if (molc%nhn(ii) == 2) then
              str = trim(str)//'R2'
            else if (molc%nhn(ii) == 3) then
              str = trim(str)//'R3'
            else if (molc%nhn(ii) == 0) then
              str = 'ammonia'
            end if
          else if (molc%hyb(ii) == 4) then
            str = trim(str)//'4+'
          end if

        case (8) !> oxygen
          if (molc%hyb(ii) == 1) then
            str = 'o'
          else if (molc%hyb(ii) == 2) then
            if (molc%nhn(ii) == 1) then
              str = trim(str)//'H'
            else if (molc%nhn(ii) == 0) then
              str = 'water'
            end if
          end if

        case (16) !> sulfur
          if (molc%hyb(ii) == 2) then
            if (molc%nhn(ii) == 1) then
              str = trim(str)//'H'
            end if
          end if

        end select
      end associate
    end do
  end subroutine atinfo_classify

!=============================================================================!
!#############################################################################!
!=============================================================================!

!> PRINTOUTS and naming

  subroutine coord_classify_print_functional(self,prch)
    implicit none
    class(coord_classify) :: self
    integer,intent(in) :: prch

    integer,allocatable :: at(:)
    integer :: ii,jj,nn

    if (.not.allocated(self%funcgroups)) return

    do ii = 1,size(self%funcgroups,1)
      nn = self%funcgroups(ii)%natms
      allocate(at(nn),source=0)
      do jj=1,nn
         at(jj) = self%at(self%funcgroups(ii)%ids(jj))
      enddo
      write (prch,'(3(1x,a))') 'functional group:', &
        & self%funcgroups(ii)%name,sumform(nn,at)
      deallocate(at)
    end do

  end subroutine coord_classify_print_functional

!=============================================================================!
!#############################################################################!
!=============================================================================!
end module molbuilder_classify_type

