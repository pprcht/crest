module test_molecular_dynamics
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use crest_testmol
  use dynamics_module
  use iomod,only:remove
  implicit none
  private

  public :: collect_mol_dynamics

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for using molecular dynamics routines in CREST
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_mol_dynamics(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
#ifdef WITH_GFNFF
    new_unittest("Compiled gfnff subproject     ",test_compiled_gfnff), &
    new_unittest("molecular dynamics (SHAKE off)",test_md_shake_off), &
    new_unittest("molecular dynamics (SHAKE on) ",test_md_shake_on), &
    new_unittest("molecular dynamics (SHAKE H)  ",test_md_shake_honly) &
#else
    new_unittest("Compiled gfnff subproject",test_compiled_gfnff,should_fail=.true.) &
#endif
    ]
!&>

  end subroutine collect_mol_dynamics

!========================================================================================!

  subroutine test_compiled_gfnff(error)
    type(error_type),allocatable,intent(out) :: error
#ifndef WITH_GFNFF
    write (*,'("       ...")') 'gfnff not compiled, expecting fail.'
    allocate (error)
#endif
  end subroutine test_compiled_gfnff

!========================================================================================!
!  The three MD tests below intentionally only set up the shared infrastructure (calcdata
!  and a test molecule) and provide placeholders for the MD-specific calls/checks.
!  Fill the marked sections with your MD driver + assertions.
!========================================================================================!

  subroutine test_md_shake_off(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    type(mddata) :: mdyn
    integer :: io
    logical :: pr

    real(wp),parameter :: e_ref = -0.6272508_wp

    !> setup calculator backend
    call sett%create('gfnff')
    call calc%add(sett)

    !> test molecule
    call get_testmol('methane',mol)

    !> MD setup
    pr = .false.
    io = 0
    call mdyn%defaults()
    mdyn%shake = .false.
    mdyn%restart = .true. !> turn on restart reading (for determinic results)
    mdyn%wrtrj = .false. !> turn off trajectory dump
    call write_fake_restart(mol,mdyn%restartfile)

    !> run
    call dynamics(mol,mdyn,calc,pr,io)

    !> cleanup
    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> checks
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,mol%energy,e_ref,thr=1e-6_wp)
    if (allocated(error)) return
  end subroutine test_md_shake_off

!========================================================================================!

  subroutine test_md_shake_on(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    type(mddata) :: mdyn
    integer :: io
    logical :: wr,pr

    real(wp),parameter :: e_ref = -0.57741556160488028_wp

    !> setup calculator backend
    call sett%create('gfnff')
    call calc%add(sett)

    !> get test molecule
    call get_testmol('methane',mol)

    !> MD setup
    pr = .false.
    io = 0
    call mdyn%defaults()
    mdyn%shake = .true.
    mdyn%restart = .true. !> turn on restart reading (for determinic results)
    mdyn%wrtrj = .false. !> turn off trajectory dump
    call write_fake_restart(mol,mdyn%restartfile)

    !> run
    call dynamics(mol,mdyn,calc,pr,io)

    !> cleanup
    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> checks
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,mol%energy,e_ref,thr=1e-6_wp)
    if (allocated(error)) return
  end subroutine test_md_shake_on

!========================================================================================!

  subroutine test_md_shake_honly(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol

    type(mddata) :: mdyn
    integer :: io
    logical :: wr,pr

    real(wp),parameter :: e_ref = -4.6456536819174667_wp

    !> setup calculator backend
    call sett%create('gfnff')
    call calc%add(sett)

    !> get test molecule
    call get_testmol('caffeine',mol)

    !> MD setup
    pr = .false.
    io = 0
    mdyn%length_ps=5.0_wp  !> shorter runtime because the mol is larger
    call mdyn%defaults()
    mdyn%shake = .true.
    mdyn%shk%shake_mode=1
    mdyn%restart = .true. !> turn on restart reading (for determinic results)_wp  !> shorter runtime because the mol is larger
    mdyn%wrtrj = .false. !> turn off trajectory dump
    call write_fake_restart(mol,mdyn%restartfile)

    !> run
    call dynamics(mol,mdyn,calc,pr,io)

    !> cleanup
    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> checks
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,mol%energy,e_ref,thr=1e-6_wp)
    if (allocated(error)) return
  end subroutine test_md_shake_honly

  subroutine write_fake_restart(mol,restartfile)
    implicit none
    type(coord),intent(in) :: mol
    character(len=:),allocatable,intent(out) :: restartfile
    integer :: ich,ii
    restartfile = 'crest_test.mdrestart'
    open (newunit=ich,file=restartfile)
    write (ich,*) 500.0_wp
    do ii = 1,mol%nat
      write (ich,'(6D22.14)') mol%xyz(1:3,ii),mol%xyz(1:3,ii)*0.0001_wp
    end do
    close (ich)
  end subroutine write_fake_restart

!========================================================================================!
!========================================================================================!
end module test_molecular_dynamics
