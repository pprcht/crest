!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2026 Christoph Plett, Philipp Pracht
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
! Ported from the xtb docking submodule (github.com/grimme-lab/xtb,
! src/prog/dock.f90), LGPL-3.0-or-later.
!================================================================================!

!> aISS docking driver, the in-process replacement for "xtb dock"
module aiss_dock
  use crest_parameters,only:wp,stdout
  use strucrd,only:coord
  use crest_calculator,only:calcdata,calculation_settings,constraint,engrad
  use aiss_param
  use aiss_iff_prepare,only:iff_precomp
  use aiss_iff_ini,only:init_iff
  use aiss_iff_energy,only:iff_e
  use aiss_search,only:aiss_docking_search
  implicit none
  private

  public :: aiss_config,aiss_wall
  public :: aiss_dock_run,aiss_add_walls

  character(len=*),parameter :: refqfile = 'aiss_refcharges'

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine aiss_dock_run(molA_in,molB_in,cfg_in,best,eint,io, &
  &                        walls,nfix,directed,ucons,iunit)
    !******************************************************************
    !* Dock molecule B onto molecule A (xtb dock A B).
    !*
    !* molA_in,molB_in : fragments in Bohr; charges/uhf come from cfg
    !* cfg_in          : run settings (copied, never modified)
    !* best            : lowest optimized complex in the CMA frame of
    !*                   the first cfg%nfrag1 atoms of A (all of A if 0)
    !* eint            : its interaction energy in kcal/mol
    !* io              : 0 on success
    !* walls           : polynomial ellipsoid walls ($wall)
    !* nfix            : atoms 1..nfix are frozen in the optimizations
    !* directed        : atom list of A for directed docking
    !* ucons           : further constraints for the complex
    !* iunit           : output unit, silent if absent
    !******************************************************************
    type(coord),intent(in)       :: molA_in,molB_in
    type(aiss_config),intent(in) :: cfg_in
    type(coord),intent(out)      :: best
    real(wp),intent(out)         :: eint
    integer,intent(out)          :: io
    type(aiss_wall),intent(in),optional :: walls(:)
    integer,intent(in),optional  :: nfix
    integer,intent(in),optional  :: directed(:)
    type(calcdata),intent(in),optional :: ucons
    integer,intent(in),optional  :: iunit

    type(aiss_config) :: cfg
    type(iff_data) :: d
    type(coord) :: molA,molB
    type(calcdata) :: calc
    real(wp) :: icoord(6),icoord0(6),r(3),e,molA_e,molB_e,etmp
    integer :: i,iu,nf
    logical :: pr

    io = 0
    eint = 0.0_wp
    iu = 0
    if (present(iunit)) iu = iunit
    pr = (iu /= 0)
    nf = 0
    if (present(nfix)) nf = nfix

    cfg = cfg_in
    if (.not.allocated(cfg%optlvl)) cfg%optlvl = 'gfn2'
    if (.not.allocated(cfg%lmolvl)) cfg%lmolvl = 'gfn2'
    if (cfg%optlvl == 'gfn0') then
      write (stdout,'(a)') 'Error: GFN0 is not supported for docking'
      io = 1
      return
    end if

    if (pr) call dock_header(iu)
    call initrand(cfg%samerand)
    call set_iff_param(iu)

    molA = molA_in
    molB = molB_in
    molA%chrg = cfg%chrg(1)
    molB%chrg = cfg%chrg(2)
    molA%uhf = cfg%uhf(1)
    molB%uhf = cfg%uhf(2)

!>--- electronic properties of both fragments
    call d%allocate(molA%nat,molB%nat)
    if (pr) write (iu,'(/,1x,a)') 'Precomputation of electronic properties'
    call iff_precomp(molA,cfg%lmolvl,cfg%etemp,d,1,etmp,io)
    if (io /= 0) return
    call iff_precomp(molB,cfg%lmolvl,cfg%etemp,d,2,etmp,io)
    if (io /= 0) return

!>--- docking CMA shift of both fragments
    if (cfg%nfrag1 .eq. 0) then
      call cmadock(molA%nat,molA%nat,molA%at,molA%xyz,r)
    else
      call cmadock(cfg%nfrag1,molA%nat,molA%at,molA%xyz,r)
    end if
    do i = 1,3
      molA%xyz(i,:) = molA%xyz(i,:)-r(i)
    end do
    call cmadock(molB%nat,molB%nat,molB%at,molB%xyz,r)
    do i = 1,3
      molB%xyz(i,:) = molB%xyz(i,:)-r(i)
    end do

!>--- IFF setup
    d%fulle = cfg%fulle
    d%debug = cfg%debug
    d%qcg = cfg%qcg
    if (present(walls)) then
      d%nwalls = size(walls)
      d%walls = walls
    end if
    icoord = 0.0_wp
    icoord0 = 0.0_wp
    call init_iff(d,cfg%nfrag1,icoord,icoord0,iu)

    if (present(directed)) then
      if (size(directed) > 0) then
        if (size(directed) > molA%nat) then
          write (stdout,'(a)') 'Error: More atoms for directed docking defined than in molecule A'
          io = 1
          return
        end if
        call setup_directed(d,directed)
      end if
    end if

!>--- cold fusion energy of the initial placement
    call iff_e(d,.false.,0,e,icoord)
    if (pr) write (iu,'(/,''Energy of cold fusion in kcal/mol:'',F8.2)') e*au

!>--- calculators, fragment references and the actual search
    call setup_complex_calc(cfg,d,molA,molB,calc,nf,walls,ucons,io)
    if (io /= 0) return
    call fragment_energies(cfg,d,molA,molB,walls,molA_e,molB_e,io)
    if (io /= 0) return

    call aiss_docking_search(cfg,d,molA,molB,calc,molA_e,molB_e,iu,best,eint,io)
    call remove_file(refqfile)
  end subroutine aiss_dock_run

!========================================================================================!

  subroutine setup_directed(d,directed)
    !*****************************************************
    !* directed docking: attractive on listed atoms, or
    !* repulsive with the distance to the list
    !*****************************************************
    type(iff_data),intent(inout) :: d
    integer,intent(in) :: directed(:)
    integer :: i,j
    real(wp) :: dist,min_dist,rep_pot

    d%ndirected = size(directed)
    d%directed_atoms = directed
    d%directed_fc = 1.0_wp
    if (allocated(d%directed_val)) deallocate (d%directed_val)
    if (allocated(d%directed_expo)) deallocate (d%directed_expo)
    allocate (d%directed_val(d%n),d%directed_expo(d%n),source=0.0_wp)

    select case (d%directed_type)
    case (p_atom_att)
      do i = 1,d%n
        if (any(i == directed)) d%directed_val(i) = attractive_pot*d%directed_fc
      end do
    case (p_atom_pot)
      do i = 1,d%n
        if (any(i == directed)) cycle
        min_dist = 0.0_wp
        do j = 1,d%ndirected
          dist = sqrt(sum((d%xyz(:,i)-d%xyz(:,directed(j)))**2))
          if (min_dist == 0.0_wp.or.dist < min_dist) min_dist = dist
        end do
        !> onset at a distance of about 4 Bohr
        rep_pot = max(0.0_wp,0.1_wp*erf(0.07_wp*min_dist-0.28_wp))
        d%directed_val(i) = rep_pot*d%directed_fc
      end do
    end select
  end subroutine setup_directed

!========================================================================================!

  subroutine setup_complex_calc(cfg,d,molA,molB,calc,nfix,walls,ucons,io)
    !*****************************************************
    !* Calculator for the complex optimizations. GFN-FF
    !* keeps A and B as separate topology fragments with
    !* fixed fragment charges (xtb: topology backup of the
    !* separated monomers).
    !*****************************************************
    type(aiss_config),intent(in) :: cfg
    type(iff_data),intent(in)    :: d
    type(coord),intent(in)       :: molA,molB
    type(calcdata),intent(out)   :: calc
    integer,intent(in)           :: nfix
    type(aiss_wall),intent(in),optional :: walls(:)
    type(calcdata),intent(in),optional :: ucons
    integer,intent(out)          :: io
    type(calculation_settings) :: clevel
    integer :: i,n
    character(len=64) :: frag1,frag2

    io = 0
    n = d%n
    call clevel%create(cfg%optlvl,chrg=molA%chrg+molB%chrg,uhf=molA%uhf+molB%uhf)
    clevel%etemp = cfg%etemp
    if (cfg%optlvl == 'gfnff') then
      write (frag1,'(i0,a,i0)') 1,'-',d%n1
      write (frag2,'(i0,a,i0)') d%n1+1,'-',n
      clevel%gff_fragments = [character(len=64) :: frag1,frag2]
      call write_refcharges(d,molA%chrg+molB%chrg)
      clevel%refcharges = refqfile
    end if
    call calc%add(clevel)
    calc%optlev = 0

    if (present(walls)) call aiss_add_walls(calc,n,walls)
    if (nfix > 0) call calc%set_freeze(n,1,nfix)
    if (present(ucons)) then
      do i = 1,ucons%nconstraints
        call calc%add(ucons%cons(i))
      end do
    end if
  end subroutine setup_complex_calc

!========================================================================================!

  subroutine fragment_energies(cfg,d,molA,molB,walls,molA_e,molB_e,io)
    !*****************************************************
    !* Unoptimized fragment energies: GFN-FF single points
    !* with the walls acting on the fragment, or the xTB
    !* precomputation energies at the same level.
    !*****************************************************
    type(aiss_config),intent(in) :: cfg
    type(iff_data),intent(in)    :: d
    type(coord),intent(in)       :: molA,molB
    type(aiss_wall),intent(in),optional :: walls(:)
    real(wp),intent(out)         :: molA_e,molB_e
    integer,intent(out)          :: io

    io = 0
    if (cfg%optlvl /= 'gfnff'.and.cfg%optlvl == cfg%lmolvl) then
      molA_e = d%etot(1)
      molB_e = d%etot(2)
      return
    end if
    call fragment_sp(molA,0,molA_e,io)
    if (io /= 0) return
    call fragment_sp(molB,molA%nat,molB_e,io)

  contains

    subroutine fragment_sp(mol,offset,energy,stat)
      type(coord),intent(in) :: mol
      integer,intent(in)     :: offset
      real(wp),intent(out)   :: energy
      integer,intent(out)    :: stat
      type(calcdata) :: fcalc
      type(calculation_settings) :: flevel
      type(coord) :: molsp
      real(wp),allocatable :: grd(:,:)
      call flevel%create(cfg%optlvl,chrg=mol%chrg,uhf=mol%uhf)
      flevel%etemp = cfg%etemp
      call fcalc%add(flevel)
      if (cfg%optlvl == 'gfnff'.and.present(walls)) call aiss_add_walls(fcalc,mol%nat,walls,offset)
      molsp = mol
      allocate (grd(3,mol%nat))
      call engrad(molsp,fcalc,energy,grd,stat)
      if (stat /= 0) write (stdout,'(a)') 'Error: fragment single point for docking failed'
    end subroutine fragment_sp

  end subroutine fragment_energies

!========================================================================================!

  subroutine aiss_add_walls(calc,nat,walls,offset)
    !*****************************************************
    !* xtb polynomial walls (|r/radius|)^30 as CREST
    !* ellipsoids with k=1, alpha=15. Wall lists use the
    !* complex numbering; offset maps them onto a fragment
    !* (atom j of the fragment is list entry j+offset).
    !*****************************************************
    type(calcdata),intent(inout) :: calc
    integer,intent(in) :: nat
    type(aiss_wall),intent(in) :: walls(:)
    integer,intent(in),optional :: offset
    type(constraint) :: cons
    logical,allocatable :: atms(:)
    integer :: i,j,k,off
    off = 0
    if (present(offset)) off = offset
    allocate (atms(nat))
    do i = 1,size(walls)
      if (allocated(walls(i)%list)) then
        atms = .false.
        do j = 1,size(walls(i)%list)
          k = walls(i)%list(j)-off
          if (k >= 1.and.k <= nat) atms(k) = .true.
        end do
      else
        atms = .true.
      end if
      if (.not.any(atms)) cycle
      call cons%ellipsoid(nat,atms,walls(i)%radius,1.0_wp,15.0_wp,.false.)
      call calc%add(cons)
    end do
  end subroutine aiss_add_walls

!========================================================================================!

  subroutine write_refcharges(d,chrg)
    !*****************************************************
    !* GFN-FF reference charges; the last atom absorbs the
    !* residue so that they sum to the total charge
    !*****************************************************
    type(iff_data),intent(in) :: d
    integer,intent(in) :: chrg
    real(wp),allocatable :: q(:)
    integer :: i,ich
    allocate (q(d%n))
    q(1:d%n1) = d%q1
    q(d%n1+1:d%n) = d%q2
    q(d%n) = q(d%n)+(real(chrg,wp)-sum(q))
    open (newunit=ich,file=refqfile)
    do i = 1,d%n
      write (ich,'(f20.12)') q(i)
    end do
    close (ich)
  end subroutine write_refcharges

!========================================================================================!

  subroutine initrand(samerand)
    logical,intent(in) :: samerand
    integer :: i,j
    integer,allocatable :: iseed(:)
    if (samerand) then
      call random_seed(size=j)
      allocate (iseed(j),source=41)
      do i = 1,j
        iseed(i) = iseed(i)+j
      end do
      call random_seed(put=iseed)
    else
      call random_seed()
    end if
  end subroutine initrand

  subroutine remove_file(fname)
    character(len=*),intent(in) :: fname
    integer :: ich,stat
    logical :: ex
    inquire (file=fname,exist=ex)
    if (.not.ex) return
    open (newunit=ich,file=fname,iostat=stat)
    if (stat == 0) close (ich,status='delete')
  end subroutine remove_file

  subroutine dock_header(iu)
    integer,intent(in) :: iu
    write (iu,'(a)') &
    & '      ------------------------------------------------------', &
    & '      |      aISS docking (xtb dock port, CREST QCG)       |', &
    & '      ------------------------------------------------------', &
    & '   Cite this work as:', &
    & '   * C. Plett, S. Grimme, Angew. Chem. Int. Ed. 2023, 62, e202214477.', &
    & '     DOI: 10.1002/anie.202214477'
  end subroutine dock_header

!========================================================================================!
end module aiss_dock
