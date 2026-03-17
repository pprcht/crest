!> symmetry_i.f90
!> Brute force symmetry analyzer - Fortran module
!>
!> Original C code: (C) 1996, 2003 S. Patchkovskii
!> Fortran conversion of the original C code
!>
!> This program is free software; you can redistribute it and/or modify
!> it under the terms of the GNU General Public License as published by
!> the Free Software Foundation; either version 2 of the License, or
!> (at your option) any later version.

! WARNING: Currently unused and untested!

module symmetry_i
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  ! Public interface
  public :: schoenflies
  public :: symmetry_element,atom_t
  public :: set_symmetry_tolerance

  !> Mathematical constants
  real(wp),parameter :: PI = 3.14159265358979323846d0
  integer,parameter :: DIMENSION = 3
  integer,parameter :: MAXPARAM = 7

  !> Atom type
  type :: atom_t
    integer :: atom_type
    real(wp) :: x(DIMENSION)
  end type atom_t

  !> Symmetry element type
  type :: symmetry_element
    integer :: transform_type  ! 1=mirror, 2=invert, 3=rotate, 4=rotate_reflect
    integer,allocatable :: transform(:)
    integer :: order
    integer :: nparam
    real(wp) :: maxdev
    real(wp) :: distance
    real(wp) :: normal(DIMENSION)
    real(wp) :: direction(DIMENSION)
  end type symmetry_element

  !> Point group type
  type :: point_group
    character(len=8) :: group_name
    character(len=64) :: symmetry_code
  end type point_group

  !> Module-level parameters (can be modified)
  real(wp),save :: ToleranceSame = 1.0d-3
  real(wp),save :: TolerancePrimary = 5.0d-2
  real(wp),save :: ToleranceFinal = 1.0d-4
  real(wp),save :: MaxOptStep = 5.0d-1
  real(wp),save :: MinOptStep = 1.0d-7
  real(wp),save :: GradientStep = 1.0d-7
  real(wp),save :: OptChangeThreshold = 1.0d-10
  integer,save :: verbose = 0
  integer,save :: MaxOptCycles = 200
  integer,save :: OptChangeHits = 5
  integer,save :: MaxAxisOrder = 20

  !> Working data
  real(wp),save :: CenterOfSomething(DIMENSION)
  real(wp),allocatable,save :: DistanceFromCenter(:)
  integer,save :: AtomsCount = 0
  type(atom_t),allocatable,save :: Atoms(:)

  !> Symmetry elements storage
  integer,save :: PlanesCount = 0
  type(symmetry_element),allocatable,save :: Planes(:)
  type(symmetry_element),allocatable,save :: MolecularPlane
  logical,save :: MolecularPlaneExists = .false.
  integer,save :: InversionCentersCount = 0
  type(symmetry_element),allocatable,save :: InversionCenters(:)
  integer,save :: NormalAxesCount = 0
  type(symmetry_element),allocatable,save :: NormalAxes(:)
  integer,save :: ImproperAxesCount = 0
  type(symmetry_element),allocatable,save :: ImproperAxes(:)
  integer,allocatable,save :: NormalAxesCounts(:)
  integer,allocatable,save :: ImproperAxesCounts(:)
  integer,save :: BadOptimization = 0
  character(len=256),save :: SymmetryCode = ""
  character(len=8),save :: MaxRotAxis = ""

  !> Statistics
  integer(8),save :: StatTotal = 0
  integer(8),save :: StatEarly = 0
  integer(8),save :: StatPairs = 0
  integer(8),save :: StatDups = 0
  integer(8),save :: StatOrder = 0
  integer(8),save :: StatOpt = 0
  integer(8),save :: StatAccept = 0

  !> Point groups table
  integer,parameter :: PointGroupsCount = 60
  type(point_group),save :: PointGroups(PointGroupsCount)
  logical,save :: PointGroupsInitialized = .false.

! ══════════════════════════════════════════════════════════════════════════════
contains    !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

  !> Initialize point groups table
  subroutine init_point_groups()
    if (PointGroupsInitialized) return

    PointGroups(1) = point_group("C1","")
    PointGroups(2) = point_group("Cs","(sigma) ")
    PointGroups(3) = point_group("Ci","(i) ")
    PointGroups(4) = point_group("C2","(C2) ")
    PointGroups(5) = point_group("C3","(C3) ")
    PointGroups(6) = point_group("C4","(C4) (C2) ")
    PointGroups(7) = point_group("C5","(C5) ")
    PointGroups(8) = point_group("C6","(C6) (C3) (C2) ")
    PointGroups(9) = point_group("C7","(C7) ")
    PointGroups(10) = point_group("C8","(C8) (C4) (C2) ")
    PointGroups(11) = point_group("D2","3*(C2) ")
    PointGroups(12) = point_group("D3","(C3) 3*(C2) ")
    PointGroups(13) = point_group("D4","(C4) 5*(C2) ")
    PointGroups(14) = point_group("D5","(C5) 5*(C2) ")
    PointGroups(15) = point_group("D6","(C6) (C3) 7*(C2) ")
    PointGroups(16) = point_group("D7","(C7) 7*(C2) ")
    PointGroups(17) = point_group("D8","(C8) (C4) 9*(C2) ")
    PointGroups(18) = point_group("C2v","(C2) 2*(sigma) ")
    PointGroups(19) = point_group("C3v","(C3) 3*(sigma) ")
    PointGroups(20) = point_group("C4v","(C4) (C2) 4*(sigma) ")
    PointGroups(21) = point_group("C5v","(C5) 5*(sigma) ")
    PointGroups(22) = point_group("C6v","(C6) (C3) (C2) 6*(sigma) ")
    PointGroups(23) = point_group("C7v","(C7) 7*(sigma) ")
    PointGroups(24) = point_group("C8v","(C8) (C4) (C2) 8*(sigma) ")
    PointGroups(25) = point_group("C2h","(i) (C2) (sigma) ")
    PointGroups(26) = point_group("C3h","(C3) (S3) (sigma) ")
    PointGroups(27) = point_group("C4h","(i) (C4) (C2) (S4) (sigma) ")
    PointGroups(28) = point_group("C5h","(C5) (S5) (sigma) ")
    PointGroups(29) = point_group("C6h","(i) (C6) (C3) (C2) (S6) (S3) (sigma) ")
    PointGroups(30) = point_group("C7h","(C7) (S7) (sigma) ")
    PointGroups(31) = point_group("C8h","(i) (C8) (C4) (C2) (S8) (S4) (sigma) ")
    PointGroups(32) = point_group("D2h","(i) 3*(C2) 3*(sigma) ")
    PointGroups(33) = point_group("D3h","(C3) 3*(C2) (S3) 4*(sigma) ")
    PointGroups(34) = point_group("D4h","(i) (C4) 5*(C2) (S4) 5*(sigma) ")
    PointGroups(35) = point_group("D5h","(C5) 5*(C2) (S5) 6*(sigma) ")
    PointGroups(36) = point_group("D6h","(i) (C6) (C3) 7*(C2) (S6) (S3) 7*(sigma) ")
    PointGroups(37) = point_group("D7h","(C7) 7*(C2) (S7) 8*(sigma) ")
    PointGroups(38) = point_group("D8h","(i) (C8) (C4) 9*(C2) (S8) (S4) 9*(sigma) ")
    PointGroups(39) = point_group("D2d","3*(C2) (S4) 2*(sigma) ")
    PointGroups(40) = point_group("D3d","(i) (C3) 3*(C2) (S6) 3*(sigma) ")
    PointGroups(41) = point_group("D4d","(C4) 5*(C2) (S8) 4*(sigma) ")
    PointGroups(42) = point_group("D5d","(i) (C5) 5*(C2) (S10) 5*(sigma) ")
    PointGroups(43) = point_group("D6d","(C6) (C3) 7*(C2) (S12) (S4) 6*(sigma) ")
    PointGroups(44) = point_group("D7d","(i) (C7) 7*(C2) (S14) 7*(sigma) ")
    PointGroups(45) = point_group("D8d","(C8) (C4) 9*(C2) (S16) 8*(sigma) ")
    PointGroups(46) = point_group("S4","(C2) (S4) ")
    PointGroups(47) = point_group("S6","(i) (C3) (S6) ")
    PointGroups(48) = point_group("S8","(C4) (C2) (S8) ")
    PointGroups(49) = point_group("T","4*(C3) 3*(C2) ")
    PointGroups(50) = point_group("Th","(i) 4*(C3) 3*(C2) 4*(S6) 3*(sigma) ")
    PointGroups(51) = point_group("Td","4*(C3) 3*(C2) 3*(S4) 6*(sigma) ")
    PointGroups(52) = point_group("O","3*(C4) 4*(C3) 9*(C2) ")
    PointGroups(53) = point_group("Oh","(i) 3*(C4) 4*(C3) 9*(C2) 4*(S6) 3*(S4) 9*(sigma) ")
    PointGroups(54) = point_group("Cinfv","(Cinf) (sigma) ")
    PointGroups(55) = point_group("Dinfh","(i) (Cinf) (C2) 2*(sigma) ")
    PointGroups(56) = point_group("I","6*(C5) 10*(C3) 15*(C2) ")
    PointGroups(57) = point_group("Ih","(i) 6*(C5) 10*(C3) 15*(C2) 6*(S10) 10*(S6) 15*(sigma) ")
    PointGroups(58) = point_group("Kh","(i) (Cinf) (sigma) ")
    PointGroups(59) = point_group("","")  ! Padding
    PointGroups(60) = point_group("","")  ! Padding

    PointGroupsInitialized = .true.
  end subroutine init_point_groups

  !> Set tolerance parameters
  subroutine set_symmetry_tolerance(tol_same,tol_primary,tol_final)
    real(wp),intent(in),optional :: tol_same,tol_primary,tol_final
    if (present(tol_same)) ToleranceSame = tol_same
    if (present(tol_primary)) TolerancePrimary = tol_primary
    if (present(tol_final)) ToleranceFinal = tol_final
  end subroutine set_symmetry_tolerance

  !> Square function
  pure real(wp) function pow2(x)
    real(wp),intent(in) :: x
    pow2 = x*x
  end function pow2

  !> Allocate a symmetry element
  subroutine alloc_symmetry_element(elem)
    type(symmetry_element),intent(out) :: elem
    integer :: i

    allocate (elem%transform(AtomsCount))
    do i = 1,AtomsCount
      elem%transform(i) = AtomsCount+1  ! Impossible value
    end do
    elem%order = 0
    elem%nparam = 0
    elem%maxdev = 0.0d0
    elem%distance = 0.0d0
    elem%normal = 0.0d0
    elem%direction = 0.0d0
    elem%transform_type = 0
  end subroutine alloc_symmetry_element

  !> Deallocate a symmetry element
  subroutine destroy_symmetry_element(elem)
    type(symmetry_element),intent(inout) :: elem
    if (allocated(elem%transform)) deallocate (elem%transform)
  end subroutine destroy_symmetry_element

  !> Mirror an atom through a plane
  subroutine mirror_atom(plane,from_atom,to_atom)
    type(symmetry_element),intent(in) :: plane
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom
    integer :: i
    real(wp) :: r

    r = plane%distance
    do i = 1,DIMENSION
      r = r-from_atom%x(i)*plane%normal(i)
    end do

    to_atom%atom_type = from_atom%atom_type
    do i = 1,DIMENSION
      to_atom%x(i) = from_atom%x(i)+2.0d0*r*plane%normal(i)
    end do
  end subroutine mirror_atom

  !> Invert an atom through a center
  subroutine invert_atom(center,from_atom,to_atom)
    type(symmetry_element),intent(in) :: center
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom
    integer :: i

    to_atom%atom_type = from_atom%atom_type
    do i = 1,DIMENSION
      to_atom%x(i) = 2.0d0*center%distance*center%normal(i)-from_atom%x(i)
    end do
  end subroutine invert_atom

  !> Rotate an atom around an axis
  subroutine rotate_atom(axis,from_atom,to_atom)
    type(symmetry_element),intent(in) :: axis
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom
    real(wp) :: x(3),y(3),a(3),b(3),c(3)
    real(wp) :: angle,a_sin,a_cos,dot_val
    integer :: i

    if (axis%order /= 0) then
      angle = 2.0d0*PI/dble(axis%order)
    else
      angle = 1.0d0
    end if
    a_sin = sin(angle)
    a_cos = cos(angle)

    do i = 1,3
      x(i) = from_atom%x(i)-axis%distance*axis%normal(i)
    end do

    dot_val = 0.0d0
    do i = 1,3
      dot_val = dot_val+x(i)*axis%direction(i)
    end do

    do i = 1,3
      a(i) = axis%direction(i)*dot_val
    end do

    do i = 1,3
      b(i) = x(i)-a(i)
    end do

    c(1) = b(2)*axis%direction(3)-b(3)*axis%direction(2)
    c(2) = b(3)*axis%direction(1)-b(1)*axis%direction(3)
    c(3) = b(1)*axis%direction(2)-b(2)*axis%direction(1)

    do i = 1,3
      y(i) = a(i)+b(i)*a_cos+c(i)*a_sin
    end do

    do i = 1,3
      to_atom%x(i) = y(i)+axis%distance*axis%normal(i)
    end do
    to_atom%atom_type = from_atom%atom_type
  end subroutine rotate_atom

  !> Rotate and reflect an atom (improper rotation)
  subroutine rotate_reflect_atom(axis,from_atom,to_atom)
    type(symmetry_element),intent(in) :: axis
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom
    real(wp) :: x(3),y(3),a(3),b(3),c(3)
    real(wp) :: angle,a_sin,a_cos,dot_val
    integer :: i

    angle = 2.0d0*PI/dble(axis%order)
    a_sin = sin(angle)
    a_cos = cos(angle)

    do i = 1,3
      x(i) = from_atom%x(i)-axis%distance*axis%normal(i)
    end do

    dot_val = 0.0d0
    do i = 1,3
      dot_val = dot_val+x(i)*axis%direction(i)
    end do

    do i = 1,3
      a(i) = axis%direction(i)*dot_val
    end do

    do i = 1,3
      b(i) = x(i)-a(i)
    end do

    c(1) = b(2)*axis%direction(3)-b(3)*axis%direction(2)
    c(2) = b(3)*axis%direction(1)-b(1)*axis%direction(3)
    c(3) = b(1)*axis%direction(2)-b(2)*axis%direction(1)

    do i = 1,3
      y(i) = -a(i)+b(i)*a_cos+c(i)*a_sin
    end do

    do i = 1,3
      to_atom%x(i) = y(i)+axis%distance*axis%normal(i)
    end do
    to_atom%atom_type = from_atom%atom_type
  end subroutine rotate_reflect_atom

  !> Transform atom based on element type
  subroutine transform_atom(elem,from_atom,to_atom)
    type(symmetry_element),intent(in) :: elem
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom

    select case (elem%transform_type)
    case (1)
      call mirror_atom(elem,from_atom,to_atom)
    case (2)
      call invert_atom(elem,from_atom,to_atom)
    case (3)
      call rotate_atom(elem,from_atom,to_atom)
    case (4)
      call rotate_reflect_atom(elem,from_atom,to_atom)
    case default
      to_atom = from_atom
    end select
  end subroutine transform_atom

  !> Establish pairs of atoms related by symmetry
  function establish_pairs(elem) result(status)
    type(symmetry_element),intent(inout) :: elem
    integer :: status
    integer :: i,j,k,best_j
    logical,allocatable :: atom_used(:)
    real(wp) :: distance,best_distance
    type(atom_t) :: symmetric

    status = 0
    allocate (atom_used(AtomsCount))
    atom_used = .false.

    do i = 1,AtomsCount
      if (elem%transform(i) > AtomsCount) then
        call transform_atom(elem,Atoms(i),symmetric)
        best_j = i
        best_distance = 2.0d0*TolerancePrimary

        do j = 1,AtomsCount
          if (Atoms(j)%atom_type /= symmetric%atom_type.or.atom_used(j)) cycle

          distance = 0.0d0
          do k = 1,DIMENSION
            distance = distance+pow2(symmetric%x(k)-Atoms(j)%x(k))
          end do
          distance = sqrt(distance)

          if (distance < best_distance) then
            best_j = j
            best_distance = distance
          end if
        end do

        if (best_distance > TolerancePrimary) then
          deallocate (atom_used)
          status = -1
          return
        end if

        elem%transform(i) = best_j
        atom_used(best_j) = .true.
      end if
    end do

    deallocate (atom_used)
  end function establish_pairs

  !> Check if transformation order is correct
  function check_transform_order(elem) result(status)
    type(symmetry_element),intent(in) :: elem
    integer :: status
    integer :: i,j,k

    status = 0

    do i = 1,AtomsCount
      if (elem%transform(i) == i) cycle

      if (elem%transform_type == 4) then  ! rotate_reflect
        j = elem%transform(i)
        if (elem%transform(j) == i) cycle
      end if

      k = elem%transform(i)
      do j = elem%order-1,1,-1
        if (k == i) then
          status = -1
          return
        end if
        k = elem%transform(k)
      end do

      if (k /= i.and.elem%transform_type == 4) then
        do j = elem%order,1,-1
          if (k == i) then
            status = -1
            return
          end if
          k = elem%transform(k)
        end do
      end if

      if (k /= i) then
        status = -1
        return
      end if
    end do
  end function check_transform_order

  !> Check if two transforms are the same
  function same_transform(a,b) result(is_same)
    type(symmetry_element),intent(in) :: a,b
    logical :: is_same
    integer :: i,j,code

    is_same = .false.

    if (a%order /= b%order.or.a%nparam /= b%nparam.or. &
        a%transform_type /= b%transform_type) return

    code = 1
    do i = 1,AtomsCount
      if (a%transform(i) /= b%transform(i)) then
        code = 0
        exit
      end if
    end do

    if (code == 0.and.a%order > 2) then
      do i = 1,AtomsCount
        j = a%transform(i)
        if (b%transform(j) /= i) return
      end do
      is_same = .true.
      return
    end if

    is_same = (code == 1)
  end function same_transform

  !> Check transform quality
  function check_transform_quality(elem) result(status)
    type(symmetry_element),intent(inout) :: elem
    integer :: status
    integer :: i,j,k
    type(atom_t) :: symmetric
    real(wp) :: r,max_r

    status = 0
    max_r = 0.0d0

    do i = 1,AtomsCount
      j = elem%transform(i)
      call transform_atom(elem,Atoms(i),symmetric)

      r = 0.0d0
      do k = 1,DIMENSION
        r = r+pow2(symmetric%x(k)-Atoms(j)%x(k))
      end do
      r = sqrt(r)

      if (r > ToleranceFinal) then
        status = -1
        return
      end if
      if (r > max_r) max_r = r
    end do

    elem%maxdev = max_r
  end function check_transform_quality

  !> Evaluate optimization target function
  function eval_optimization_target_function(elem,finish) result(target)
    type(symmetry_element),intent(inout) :: elem
    logical,intent(out),optional :: finish
    real(wp) :: target
    integer :: i,j,k
    type(atom_t) :: symmetric
    real(wp) :: r,maxr

    ! Normalize normal vector
    if (elem%nparam >= 4) then
      r = 0.0d0
      do k = 1,DIMENSION
        r = r+elem%normal(k)*elem%normal(k)
      end do
      r = sqrt(r)
      if (r < ToleranceSame) then
        write (*,*) "Normal collapsed!"
        stop
      end if
      elem%normal = elem%normal/r
      if (elem%distance < 0.0d0) then
        elem%distance = -elem%distance
        elem%normal = -elem%normal
      end if
    end if

    ! Normalize direction vector
    if (elem%nparam >= 7) then
      r = 0.0d0
      do k = 1,DIMENSION
        r = r+elem%direction(k)*elem%direction(k)
      end do
      r = sqrt(r)
      if (r < ToleranceSame) then
        write (*,*) "Direction collapsed!"
        stop
      end if
      elem%direction = elem%direction/r
    end if

    target = 0.0d0
    maxr = 0.0d0

    do i = 1,AtomsCount
      call transform_atom(elem,Atoms(i),symmetric)
      j = elem%transform(i)

      r = 0.0d0
      do k = 1,DIMENSION
        r = r+pow2(Atoms(j)%x(k)-symmetric%x(k))
      end do
      if (r > maxr) maxr = r
      target = target+r
    end do

    if (present(finish)) then
      finish = (sqrt(maxr) < ToleranceFinal)
    end if
  end function eval_optimization_target_function

  !> Get parameters from element
  subroutine get_params(elem,values)
    type(symmetry_element),intent(in) :: elem
    real(wp),intent(out) :: values(MAXPARAM)

    values(1) = elem%distance
    values(2:4) = elem%normal(1:3)
    if (elem%nparam >= 7) then
      values(5:7) = elem%direction(1:3)
    end if
  end subroutine get_params

  !> Set parameters to element
  subroutine set_params(elem,values)
    type(symmetry_element),intent(inout) :: elem
    real(wp),intent(in) :: values(MAXPARAM)

    elem%distance = values(1)
    elem%normal(1:3) = values(2:4)
    if (elem%nparam >= 7) then
      elem%direction(1:3) = values(5:7)
    end if
  end subroutine set_params

  !> Optimize transformation parameters
  subroutine optimize_transformation_params(elem)
    type(symmetry_element),intent(inout) :: elem
    real(wp) :: values(MAXPARAM),grad(MAXPARAM),force(MAXPARAM),step(MAXPARAM)
    real(wp) :: f,fold,fnew,fnew2,fdn,fup,snorm
    real(wp) :: a,b,x
    integer :: vars,cycle,i,hits
    logical :: finish

    vars = elem%nparam
    if (vars > MAXPARAM) then
      write (*,*) "Catastrophe in optimize_transformation_params!"
      stop
    end if

    f = 0.0d0
    cycle = 0
    hits = 0

    do
      fold = f
      f = eval_optimization_target_function(elem,finish)

      if (finish) exit

      if (cycle > 0) then
        if (abs(f-fold) > OptChangeThreshold) then
          hits = 0
        else
          hits = hits+1
        end if
        if (hits >= OptChangeHits) exit
      end if

      call get_params(elem,values)

      ! Calculate gradient and force constants
      do i = 1,vars
        values(i) = values(i)-GradientStep
        call set_params(elem,values)
        fdn = eval_optimization_target_function(elem)

        values(i) = values(i)+2.0d0*GradientStep
        call set_params(elem,values)
        fup = eval_optimization_target_function(elem)

        values(i) = values(i)-GradientStep
        grad(i) = (fup-fdn)/(2.0d0*GradientStep)
        force(i) = (fup+fdn-2.0d0*f)/(GradientStep*GradientStep)
      end do

      ! Quasi-Newton step
      snorm = 0.0d0
      do i = 1,vars
        if (force(i) < 0.0d0) force(i) = -force(i)
        if (force(i) < 1.0d-3) force(i) = 1.0d-3
        if (force(i) > 1.0d3) force(i) = 1.0d3
        step(i) = -grad(i)/force(i)
        snorm = snorm+step(i)*step(i)
      end do
      snorm = sqrt(snorm)

      if (snorm > MaxOptStep) then
        step = step*MaxOptStep/snorm
        snorm = MaxOptStep
      end if

      do while (snorm > MinOptStep)
        values = values+step
        call set_params(elem,values)
        fnew = eval_optimization_target_function(elem)

        if (fnew < f) exit

        values = values-step
        step = step/2.0d0
        call set_params(elem,values)
        snorm = snorm/2.0d0
      end do

      ! Quadratic interpolation
      if (snorm > MinOptStep.and.snorm < MaxOptStep/2.0d0) then
        values = values+step
        call set_params(elem,values)
        fnew2 = eval_optimization_target_function(elem)
        values = values-2.0d0*step

        a = (4.0d0*f-fnew2-3.0d0*fnew)/2.0d0
        b = (f+fnew2-2.0d0*fnew)/2.0d0

        if (b > 0.0d0) then
          x = -a/(2.0d0*b)
          if (x > 0.2d0.and.x < 1.8d0) then
            values = values+x*step
          else
            b = 0.0d0
          end if
        end if

        if (b <= 0.0d0) then
          if (fnew2 < fnew) then
            values = values+2.0d0*step
          else
            values = values+step
          end if
        end if
        call set_params(elem,values)
      end if

      cycle = cycle+1
      if (snorm <= MinOptStep.or.cycle >= MaxOptCycles) exit
    end do

    f = eval_optimization_target_function(elem)
    if (cycle >= MaxOptCycles) BadOptimization = 1
  end subroutine optimize_transformation_params

  !> Refine symmetry element
  function refine_symmetry_element(elem,build_table) result(status)
    type(symmetry_element),intent(inout) :: elem
    logical,intent(in) :: build_table
    integer :: status
    integer :: i

    status = 0

    if (build_table) then
      if (establish_pairs(elem) < 0) then
        StatPairs = StatPairs+1
        status = -1
        return
      end if
    end if

    ! Check for duplicates
    do i = 1,PlanesCount
      if (same_transform(Planes(i),elem)) then
        StatDups = StatDups+1
        status = -1
        return
      end if
    end do

    do i = 1,InversionCentersCount
      if (same_transform(InversionCenters(i),elem)) then
        StatDups = StatDups+1
        status = -1
        return
      end if
    end do

    do i = 1,NormalAxesCount
      if (same_transform(NormalAxes(i),elem)) then
        StatDups = StatDups+1
        status = -1
        return
      end if
    end do

    do i = 1,ImproperAxesCount
      if (same_transform(ImproperAxes(i),elem)) then
        StatDups = StatDups+1
        status = -1
        return
      end if
    end do

    if (check_transform_order(elem) < 0) then
      StatOrder = StatOrder+1
      status = -1
      return
    end if

    call optimize_transformation_params(elem)

    if (check_transform_quality(elem) < 0) then
      StatOpt = StatOpt+1
      status = -1
      return
    end if

    StatAccept = StatAccept+1
  end function refine_symmetry_element

  !> Initialize mirror plane
  subroutine init_mirror_plane(i,j,plane,success)
    integer,intent(in) :: i,j
    type(symmetry_element),intent(out) :: plane
    logical,intent(out) :: success
    real(wp) :: dx(DIMENSION),midpoint(DIMENSION),rab,r
    integer :: k

    success = .false.
    StatTotal = StatTotal+1

    call alloc_symmetry_element(plane)
    plane%transform_type = 1  ! mirror
    plane%order = 2
    plane%nparam = 4

    rab = 0.0d0
    do k = 1,DIMENSION
      dx(k) = Atoms(i)%x(k)-Atoms(j)%x(k)
      midpoint(k) = (Atoms(i)%x(k)+Atoms(j)%x(k))/2.0d0
      rab = rab+dx(k)*dx(k)
    end do
    rab = sqrt(rab)

    if (rab < ToleranceSame) then
      call destroy_symmetry_element(plane)
      return
    end if

    r = 0.0d0
    do k = 1,DIMENSION
      plane%normal(k) = dx(k)/rab
      r = r+midpoint(k)*plane%normal(k)
    end do

    if (r < 0.0d0) then
      r = -r
      plane%normal = -plane%normal
    end if
    plane%distance = r

    if (refine_symmetry_element(plane,.true.) < 0) then
      call destroy_symmetry_element(plane)
      return
    end if

    success = .true.
  end subroutine init_mirror_plane

  !> Initialize ultimate (whole-molecule) plane
  subroutine init_ultimate_plane(plane,success)
    type(symmetry_element),intent(out) :: plane
    logical,intent(out) :: success
    real(wp) :: d0(DIMENSION),d1(DIMENSION),d2(DIMENSION),p(DIMENSION)
    real(wp) :: r,s0,s1,s2
    real(wp),pointer :: d(:)
    integer :: i,j,k

    success = .false.
    StatTotal = StatTotal+1

    call alloc_symmetry_element(plane)
    plane%transform_type = 1
    plane%order = 1
    plane%nparam = 4

    d0 = 0.0d0; d1 = 0.0d0; d2 = 0.0d0
    d0(1) = 1.0d0; d1(2) = 1.0d0; d2(3) = 1.0d0

    do i = 2,AtomsCount
      do j = 1,i-1
        r = 0.0d0
        do k = 1,DIMENSION
          p(k) = Atoms(i)%x(k)-Atoms(j)%x(k)
          r = r+p(k)*p(k)
        end do
        r = sqrt(r)

        s0 = 0.0d0; s1 = 0.0d0; s2 = 0.0d0
        do k = 1,DIMENSION
          p(k) = p(k)/r
          s0 = s0+p(k)*d0(k)
          s1 = s1+p(k)*d1(k)
          s2 = s2+p(k)*d2(k)
        end do

        do k = 1,DIMENSION
          d0(k) = d0(k)-s0*p(k)
          d1(k) = d1(k)-s1*p(k)
          d2(k) = d2(k)-s2*p(k)
        end do
      end do
    end do

    s0 = sum(d0)
    s1 = sum(d1)
    s2 = sum(d2)

    if (s0 >= s1.and.s0 >= s2) then
      plane%normal = d0
    else if (s1 >= s0.and.s1 >= s2) then
      plane%normal = d1
    else
      plane%normal = d2
    end if

    r = sqrt(sum(plane%normal**2))
    if (r > 0.0d0) then
      plane%normal = plane%normal/r
    else
      plane%normal = [1.0d0,0.0d0,0.0d0]
    end if

    r = dot_product(CenterOfSomething,plane%normal)
    plane%distance = r

    do k = 1,AtomsCount
      plane%transform(k) = k
    end do

    if (refine_symmetry_element(plane,.false.) < 0) then
      call destroy_symmetry_element(plane)
      return
    end if

    success = .true.
  end subroutine init_ultimate_plane

  !> Initialize inversion center
  subroutine init_inversion_center(center,success)
    type(symmetry_element),intent(out) :: center
    logical,intent(out) :: success
    real(wp) :: r
    integer :: k

    success = .false.
    StatTotal = StatTotal+1

    call alloc_symmetry_element(center)
    center%transform_type = 2  ! invert
    center%order = 2
    center%nparam = 4

    r = sqrt(sum(CenterOfSomething**2))

    if (r > 0.0d0) then
      center%normal = CenterOfSomething/r
    else
      center%normal = [1.0d0,0.0d0,0.0d0]
    end if
    center%distance = r

    if (refine_symmetry_element(center,.true.) < 0) then
      call destroy_symmetry_element(center)
      return
    end if

    success = .true.
  end subroutine init_inversion_center

  !> Initialize ultimate (infinity) axis
  subroutine init_ultimate_axis(axis,success)
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: dir(DIMENSION),rel(DIMENSION),s
    integer :: i,k

    success = .false.
    StatTotal = StatTotal+1

    call alloc_symmetry_element(axis)
    axis%transform_type = 3  ! rotate
    axis%order = 0
    axis%nparam = 7

    dir = 0.0d0
    do i = 1,AtomsCount
      s = 0.0d0
      do k = 1,DIMENSION
        rel(k) = Atoms(i)%x(k)-CenterOfSomething(k)
        s = s+rel(k)*dir(k)
      end do
      if (s >= 0.0d0) then
        dir = dir+rel
      else
        dir = dir-rel
      end if
    end do

    s = sqrt(sum(dir**2))
    if (s > 0.0d0) then
      axis%direction = dir/s
    else
      axis%direction = [1.0d0,0.0d0,0.0d0]
    end if

    s = sqrt(sum(CenterOfSomething**2))
    if (s > 0.0d0) then
      axis%normal = CenterOfSomething/s
    else
      axis%normal = [1.0d0,0.0d0,0.0d0]
    end if
    axis%distance = s

    do k = 1,AtomsCount
      axis%transform(k) = k
    end do

    if (refine_symmetry_element(axis,.false.) < 0) then
      call destroy_symmetry_element(axis)
      return
    end if

    success = .true.
  end subroutine init_ultimate_axis

  !> Initialize C2 axis
  subroutine init_c2_axis(i,j,support,axis,success)
    integer,intent(in) :: i,j
    real(wp),intent(in) :: support(DIMENSION)
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: ris,rjs,r,center(DIMENSION)
    integer :: k

    success = .false.
    StatTotal = StatTotal+1

    ! Quick sanity check
    ris = 0.0d0
    rjs = 0.0d0
    do k = 1,DIMENSION
      ris = ris+pow2(Atoms(i)%x(k)-support(k))
      rjs = rjs+pow2(Atoms(j)%x(k)-support(k))
    end do
    ris = sqrt(ris)
    rjs = sqrt(rjs)

    if (abs(ris-rjs) > TolerancePrimary) then
      StatEarly = StatEarly+1
      return
    end if

    call alloc_symmetry_element(axis)
    axis%transform_type = 3  ! rotate
    axis%order = 2
    axis%nparam = 7

    r = sqrt(sum(CenterOfSomething**2))
    if (r > 0.0d0) then
      axis%normal = CenterOfSomething/r
    else
      axis%normal = [1.0d0,0.0d0,0.0d0]
    end if
    axis%distance = r

    r = 0.0d0
    do k = 1,DIMENSION
      center(k) = (Atoms(i)%x(k)+Atoms(j)%x(k))/2.0d0-support(k)
      r = r+center(k)*center(k)
    end do
    r = sqrt(r)

    if (r <= TolerancePrimary) then
      ! C2 is underdefined
      if (MolecularPlaneExists) then
        axis%direction = MolecularPlane%normal
      else
        do k = 1,DIMENSION
          center(k) = Atoms(i)%x(k)-Atoms(j)%x(k)
        end do
        if (abs(center(3))+abs(center(2)) > ToleranceSame) then
          axis%direction = [0.0d0,center(3),-center(2)]
        else
          axis%direction = [-center(3),0.0d0,center(1)]
        end if
        r = sqrt(sum(axis%direction**2))
        axis%direction = axis%direction/r
      end if
    else
      axis%direction = center/r
    end if

    if (refine_symmetry_element(axis,.true.) < 0) then
      call destroy_symmetry_element(axis)
      return
    end if

    success = .true.
  end subroutine init_c2_axis

  !> Initialize axis parameters from three points
  subroutine init_axis_parameters(a,b,c,axis,success)
    real(wp),intent(in) :: a(3),b(3),c(3)
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: ra,rb,rc,rab,rbc,rac,r,angle
    integer :: i,order,sign_val

    success = .false.

    ra = sqrt(sum(a**2))
    rb = sqrt(sum(b**2))
    rc = sqrt(sum(c**2))

    if (abs(ra-rb) > TolerancePrimary.or. &
        abs(ra-rc) > TolerancePrimary.or. &
        abs(rb-rc) > TolerancePrimary) then
      StatEarly = StatEarly+1
      return
    end if

    rab = sqrt(sum((a-b)**2))
    rac = sqrt(sum((a-c)**2))
    rbc = sqrt(sum((c-b)**2))

    if (abs(rab-rbc) > TolerancePrimary) then
      StatEarly = StatEarly+1
      return
    end if

    if (rab <= ToleranceSame.or.rbc <= ToleranceSame.or.rac <= ToleranceSame) then
      StatEarly = StatEarly+1
      return
    end if

    rab = (rab+rbc)/2.0d0
    angle = PI-2.0d0*asin(rac/(2.0d0*rab))

    if (abs(angle) <= PI/(MaxAxisOrder+1)) then
      StatEarly = StatEarly+1
      return
    end if

    order = nint((2.0d0*PI)/angle)
    if (order <= 2.or.order > MaxAxisOrder) then
      StatEarly = StatEarly+1
      return
    end if

    call alloc_symmetry_element(axis)
    axis%order = order
    axis%nparam = 7

    r = sqrt(sum(CenterOfSomething**2))
    if (r > 0.0d0) then
      axis%normal = CenterOfSomething/r
    else
      axis%normal = [1.0d0,0.0d0,0.0d0]
    end if
    axis%distance = r

    ! Cross product for direction
    axis%direction(1) = (b(2)-a(2))*(c(3)-b(3))-(b(3)-a(3))*(c(2)-b(2))
    axis%direction(2) = (b(3)-a(3))*(c(1)-b(1))-(b(1)-a(1))*(c(3)-b(3))
    axis%direction(3) = (b(1)-a(1))*(c(2)-b(2))-(b(2)-a(2))*(c(1)-b(1))

    ! Select direction so first non-zero component is positive
    sign_val = 0
    if (axis%direction(1) < 0.0d0) then
      sign_val = 1
    else if (axis%direction(1) == 0.0d0) then
      if (axis%direction(2) < 0.0d0) then
        sign_val = 1
      else if (axis%direction(2) == 0.0d0) then
        if (axis%direction(3) < 0.0d0) sign_val = 1
      end if
    end if

    if (sign_val == 1) axis%direction = -axis%direction

    r = sqrt(sum(axis%direction**2))
    axis%direction = axis%direction/r

    success = .true.
  end subroutine init_axis_parameters

  !> Initialize higher-order axis
  subroutine init_higher_axis(ia,ib,ic,axis,success)
    integer,intent(in) :: ia,ib,ic
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: a(DIMENSION),b(DIMENSION),c(DIMENSION)
    integer :: i

    success = .false.
    StatTotal = StatTotal+1

    do i = 1,DIMENSION
      a(i) = Atoms(ia)%x(i)-CenterOfSomething(i)
      b(i) = Atoms(ib)%x(i)-CenterOfSomething(i)
      c(i) = Atoms(ic)%x(i)-CenterOfSomething(i)
    end do

    call init_axis_parameters(a,b,c,axis,success)
    if (.not.success) return

    axis%transform_type = 3  ! rotate

    if (refine_symmetry_element(axis,.true.) < 0) then
      call destroy_symmetry_element(axis)
      success = .false.
      return
    end if

    success = .true.
  end subroutine init_higher_axis

  !> Initialize improper axis
  subroutine init_improper_axis(ia,ib,ic,axis,success)
    integer,intent(in) :: ia,ib,ic
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: a(DIMENSION),b(DIMENSION),c(DIMENSION)
    real(wp) :: centerpoint(DIMENSION),r
    integer :: i

    success = .false.
    StatTotal = StatTotal+1

    do i = 1,DIMENSION
      a(i) = Atoms(ia)%x(i)-CenterOfSomething(i)
      b(i) = Atoms(ib)%x(i)-CenterOfSomething(i)
      c(i) = Atoms(ic)%x(i)-CenterOfSomething(i)
    end do

    r = 0.0d0
    do i = 1,DIMENSION
      centerpoint(i) = a(i)+c(i)+2.0d0*b(i)
      r = r+centerpoint(i)*centerpoint(i)
    end do
    r = sqrt(r)

    if (r <= ToleranceSame) then
      StatEarly = StatEarly+1
      return
    end if

    centerpoint = centerpoint/r
    r = dot_product(centerpoint,b)
    b = 2.0d0*r*centerpoint-b

    call init_axis_parameters(a,b,c,axis,success)
    if (.not.success) return

    axis%transform_type = 4  ! rotate_reflect

    if (refine_symmetry_element(axis,.true.) < 0) then
      call destroy_symmetry_element(axis)
      success = .false.
      return
    end if

    success = .true.
  end subroutine init_improper_axis

  !> Find center of something (centroid)
  subroutine find_center_of_something()
    integer :: i,j
    real(wp) :: coord_sum(DIMENSION),r

    coord_sum = 0.0d0
    do i = 1,AtomsCount
      coord_sum = coord_sum+Atoms(i)%x
    end do
    CenterOfSomething = coord_sum/dble(AtomsCount)

    if (allocated(DistanceFromCenter)) deallocate (DistanceFromCenter)
    allocate (DistanceFromCenter(AtomsCount))

    do i = 1,AtomsCount
      r = 0.0d0
      do j = 1,DIMENSION
        r = r+pow2(Atoms(i)%x(j)-CenterOfSomething(j))
      end do
      DistanceFromCenter(i) = r
    end do
  end subroutine find_center_of_something

  !> Add plane to planes array
  subroutine add_plane(plane)
    type(symmetry_element),intent(in) :: plane
    type(symmetry_element),allocatable :: temp(:)

    PlanesCount = PlanesCount+1
    if (allocated(Planes)) then
      allocate (temp(PlanesCount))
      temp(1:PlanesCount-1) = Planes
      temp(PlanesCount) = plane
      call move_alloc(temp,Planes)
    else
      allocate (Planes(1))
      Planes(1) = plane
    end if
  end subroutine add_plane

  !> Add normal axis to array
  subroutine add_normal_axis(axis)
    type(symmetry_element),intent(in) :: axis
    type(symmetry_element),allocatable :: temp(:)

    NormalAxesCount = NormalAxesCount+1
    if (allocated(NormalAxes)) then
      allocate (temp(NormalAxesCount))
      temp(1:NormalAxesCount-1) = NormalAxes
      temp(NormalAxesCount) = axis
      call move_alloc(temp,NormalAxes)
    else
      allocate (NormalAxes(1))
      NormalAxes(1) = axis
    end if
  end subroutine add_normal_axis

  !> Add improper axis to array
  subroutine add_improper_axis(axis)
    type(symmetry_element),intent(in) :: axis
    type(symmetry_element),allocatable :: temp(:)

    ImproperAxesCount = ImproperAxesCount+1
    if (allocated(ImproperAxes)) then
      allocate (temp(ImproperAxesCount))
      temp(1:ImproperAxesCount-1) = ImproperAxes
      temp(ImproperAxesCount) = axis
      call move_alloc(temp,ImproperAxes)
    else
      allocate (ImproperAxes(1))
      ImproperAxes(1) = axis
    end if
  end subroutine add_improper_axis

  !> Find planes of symmetry
  subroutine find_planes()
    integer :: i,j
    type(symmetry_element) :: plane
    logical :: success

    call init_ultimate_plane(plane,success)
    if (success) then
      if (.not.allocated(MolecularPlane)) allocate (MolecularPlane)
      MolecularPlane = plane
      MolecularPlaneExists = .true.
      call add_plane(plane)
    end if

    do i = 2,AtomsCount
      do j = 1,i-1
        if (Atoms(i)%atom_type /= Atoms(j)%atom_type) cycle

        call init_mirror_plane(i,j,plane,success)
        if (success) call add_plane(plane)
      end do
    end do
  end subroutine find_planes

  !> Find inversion centers
  subroutine find_inversion_centers()
    type(symmetry_element) :: center
    logical :: success

    call init_inversion_center(center,success)
    if (success) then
      InversionCentersCount = 1
      allocate (InversionCenters(1))
      InversionCenters(1) = center
    end if
  end subroutine find_inversion_centers

  !> Find infinity axis
  subroutine find_infinity_axis()
    type(symmetry_element) :: axis
    logical :: success

    call init_ultimate_axis(axis,success)
    if (success) call add_normal_axis(axis)
  end subroutine find_infinity_axis

  !> Find C2 axes
  subroutine find_c2_axes()
    integer :: i,j,k,l,m
    real(wp) :: center(DIMENSION),r
    real(wp),allocatable :: distances(:)
    type(symmetry_element) :: axis
    logical :: success

    allocate (distances(AtomsCount))

    do i = 2,AtomsCount
      do j = 1,i-1
        if (Atoms(i)%atom_type /= Atoms(j)%atom_type) cycle
        if (abs(DistanceFromCenter(i)-DistanceFromCenter(j)) > TolerancePrimary) cycle

        ! Try using CenterOfSomething
        r = 0.0d0
        do k = 1,DIMENSION
          center(k) = (Atoms(i)%x(k)+Atoms(j)%x(k))/2.0d0
          r = r+pow2(center(k)-CenterOfSomething(k))
        end do
        r = sqrt(r)

        if (r > 5.0d0*TolerancePrimary) then
          call init_c2_axis(i,j,CenterOfSomething,axis,success)
          if (success) call add_normal_axis(axis)
          cycle
        end if

        ! Try through atoms
        do k = 1,AtomsCount
          call init_c2_axis(i,j,Atoms(k)%x,axis,success)
          if (success) call add_normal_axis(axis)
        end do

        ! Calculate distances for prescreening
        do k = 1,AtomsCount
          r = 0.0d0
          do l = 1,DIMENSION
            r = r+pow2(Atoms(k)%x(l)-center(l))
          end do
          distances(k) = sqrt(r)
        end do

        ! Try through midpoints of atom pairs
        do k = 1,AtomsCount
          do l = 1,AtomsCount
            if (Atoms(k)%atom_type /= Atoms(l)%atom_type) cycle
            if (abs(DistanceFromCenter(k)-DistanceFromCenter(l)) > TolerancePrimary.or. &
                abs(distances(k)-distances(l)) > TolerancePrimary) cycle

            do m = 1,DIMENSION
              center(m) = (Atoms(k)%x(m)+Atoms(l)%x(m))/2.0d0
            end do

            call init_c2_axis(i,j,center,axis,success)
            if (success) call add_normal_axis(axis)
          end do
        end do
      end do
    end do

    deallocate (distances)
  end subroutine find_c2_axes

  !> Find higher-order axes
  subroutine find_higher_axes()
    integer :: i,j,k
    type(symmetry_element) :: axis
    logical :: success

    do i = 1,AtomsCount
      do j = i+1,AtomsCount
        if (Atoms(i)%atom_type /= Atoms(j)%atom_type) cycle
        if (abs(DistanceFromCenter(i)-DistanceFromCenter(j)) > TolerancePrimary) cycle

        do k = 1,AtomsCount
          if (Atoms(i)%atom_type /= Atoms(k)%atom_type) cycle
          if (abs(DistanceFromCenter(i)-DistanceFromCenter(k)) > TolerancePrimary.or. &
              abs(DistanceFromCenter(j)-DistanceFromCenter(k)) > TolerancePrimary) cycle

          call init_higher_axis(i,j,k,axis,success)
          if (success) call add_normal_axis(axis)
        end do
      end do
    end do
  end subroutine find_higher_axes

  !> Find improper axes
  subroutine find_improper_axes()
    integer :: i,j,k
    type(symmetry_element) :: axis
    logical :: success

    do i = 1,AtomsCount
      do j = i+1,AtomsCount
        do k = 1,AtomsCount
          call init_improper_axis(i,j,k,axis,success)
          if (success) call add_improper_axis(axis)
        end do
      end do
    end do
  end subroutine find_improper_axes

  !> Find all symmetry elements
  subroutine find_symmetry_elements()
    call find_center_of_something()
    call find_inversion_centers()
    call find_planes()
    call find_infinity_axis()
    call find_c2_axes()
    call find_higher_axes()
    call find_improper_axes()
  end subroutine find_symmetry_elements

  !> Compare axes for sorting
  function compare_axes(a,b) result(cmp)
    type(symmetry_element),intent(in) :: a,b
    integer :: cmp
    integer :: order_a,order_b

    order_a = a%order
    order_b = b%order
    if (order_a == 0) order_a = 10000
    if (order_b == 0) order_b = 10000

    cmp = order_b-order_a
    if (cmp /= 0) return

    if (a%maxdev > b%maxdev) then
      cmp = -1
    else if (a%maxdev < b%maxdev) then
      cmp = 1
    else
      cmp = 0
    end if
  end function compare_axes

  !> Sort symmetry elements (simple bubble sort)
  subroutine sort_symmetry_elements()
    integer :: i,j
    type(symmetry_element) :: temp

    ! Sort planes
    do i = 1,PlanesCount-1
      do j = i+1,PlanesCount
        if (compare_axes(Planes(i),Planes(j)) < 0) then
          temp = Planes(i)
          Planes(i) = Planes(j)
          Planes(j) = temp
        end if
      end do
    end do

    ! Sort normal axes
    do i = 1,NormalAxesCount-1
      do j = i+1,NormalAxesCount
        if (compare_axes(NormalAxes(i),NormalAxes(j)) < 0) then
          temp = NormalAxes(i)
          NormalAxes(i) = NormalAxes(j)
          NormalAxes(j) = temp
        end if
      end do
    end do

    ! Sort improper axes
    do i = 1,ImproperAxesCount-1
      do j = i+1,ImproperAxesCount
        if (compare_axes(ImproperAxes(i),ImproperAxes(j)) < 0) then
          temp = ImproperAxes(i)
          ImproperAxes(i) = ImproperAxes(j)
          ImproperAxes(j) = temp
        end if
      end do
    end do
  end subroutine sort_symmetry_elements

  !> Summarize symmetry elements
  subroutine summarize_symmetry_elements()
    integer :: i

    if (allocated(NormalAxesCounts)) deallocate (NormalAxesCounts)
    if (allocated(ImproperAxesCounts)) deallocate (ImproperAxesCounts)

    allocate (NormalAxesCounts(0:MaxAxisOrder))
    allocate (ImproperAxesCounts(0:MaxAxisOrder))

    NormalAxesCounts = 0
    ImproperAxesCounts = 0

    do i = 1,NormalAxesCount
      NormalAxesCounts(NormalAxes(i)%order) = NormalAxesCounts(NormalAxes(i)%order)+1
    end do

    do i = 1,ImproperAxesCount
      ImproperAxesCounts(ImproperAxes(i)%order) = ImproperAxesCounts(ImproperAxes(i)%order)+1
    end do
  end subroutine summarize_symmetry_elements

  !> Report symmetry elements brief
  subroutine report_symmetry_elements_brief()
    integer :: i
    character(len=32) :: buf

    SymmetryCode = ""

    if (PlanesCount+NormalAxesCount+ImproperAxesCount+InversionCentersCount > 0) then
      if (InversionCentersCount > 0) SymmetryCode = trim(SymmetryCode)//"(i) "

      if (NormalAxesCounts(0) == 1) then
        SymmetryCode = trim(SymmetryCode)//"(Cinf) "
      else if (NormalAxesCounts(0) > 1) then
        write (buf,'(I0,A)') NormalAxesCounts(0),"*(Cinf) "
        SymmetryCode = trim(SymmetryCode)//trim(buf)
      end if

      do i = MaxAxisOrder,2,-1
        if (NormalAxesCounts(i) == 1) then
          write (buf,'(A,I0,A)') "(C",i,") "
          SymmetryCode = trim(SymmetryCode)//trim(buf)
        else if (NormalAxesCounts(i) > 1) then
          write (buf,'(I0,A,I0,A)') NormalAxesCounts(i),"*(C",i,") "
          SymmetryCode = trim(SymmetryCode)//trim(buf)
        end if
      end do

      do i = MaxAxisOrder,2,-1
        if (ImproperAxesCounts(i) == 1) then
          write (buf,'(A,I0,A)') "(S",i,") "
          SymmetryCode = trim(SymmetryCode)//trim(buf)
        else if (ImproperAxesCounts(i) > 1) then
          write (buf,'(I0,A,I0,A)') ImproperAxesCounts(i),"*(S",i,") "
          SymmetryCode = trim(SymmetryCode)//trim(buf)
        end if
      end do

      if (PlanesCount == 1) then
        SymmetryCode = trim(SymmetryCode)//"(sigma) "
      else if (PlanesCount > 1) then
        write (buf,'(I0,A)') PlanesCount,"*(sigma) "
        SymmetryCode = trim(SymmetryCode)//trim(buf)
      end if
    end if
  end subroutine report_symmetry_elements_brief

  !> Report highest rotation axis only
  subroutine report_symmetry_elements_brief_conly()
    integer :: i
    character(len=8) :: buf

    MaxRotAxis = ""

    if (PlanesCount+NormalAxesCount+ImproperAxesCount+InversionCentersCount > 0) then
      do i = MaxAxisOrder,2,-1
        if (NormalAxesCounts(i) >= 1) then
          write (buf,'(A,I0)') "C",i
          MaxRotAxis = trim(buf)
          return
        end if
      end do
    end if
  end subroutine report_symmetry_elements_brief_conly

  !> Identify point group
  function identify_point_group() result(last_matching)
    integer :: last_matching
    integer :: i,matching_count

    call init_point_groups()

    last_matching = -1
    matching_count = 0

    do i = 1,PointGroupsCount
      if (len_trim(PointGroups(i)%group_name) == 0) cycle
      if (trim(SymmetryCode) == trim(PointGroups(i)%symmetry_code)) then
        last_matching = i
        matching_count = matching_count+1
      end if
    end do

    if (matching_count == 0) then
      last_matching = -1
    else if (matching_count > 1) then
      last_matching = -1
    end if
  end function identify_point_group

  !> Reset module state
  subroutine reset_state()
    PlanesCount = 0
    InversionCentersCount = 0
    NormalAxesCount = 0
    ImproperAxesCount = 0
    BadOptimization = 0
    SymmetryCode = ""
    MaxRotAxis = ""
    MolecularPlaneExists = .false.

    StatTotal = 0
    StatEarly = 0
    StatPairs = 0
    StatDups = 0
    StatOrder = 0
    StatOpt = 0
    StatAccept = 0

    if (allocated(Planes)) deallocate (Planes)
    if (allocated(MolecularPlane)) deallocate (MolecularPlane)
    if (allocated(InversionCenters)) deallocate (InversionCenters)
    if (allocated(NormalAxes)) deallocate (NormalAxes)
    if (allocated(ImproperAxes)) deallocate (ImproperAxes)
    if (allocated(NormalAxesCounts)) deallocate (NormalAxesCounts)
    if (allocated(ImproperAxesCounts)) deallocate (ImproperAxesCounts)
    if (allocated(DistanceFromCenter)) deallocate (DistanceFromCenter)
    if (allocated(Atoms)) deallocate (Atoms)
  end subroutine reset_state

  !> Main entry point: determine Schoenflies symbol
  subroutine schoenflies(natoms,attype,coord,symbol,paramar)
    integer,intent(in) :: natoms
    integer,intent(in) :: attype(natoms)
    real(wp),intent(in) :: coord(3,natoms)
    character(len=*),intent(out) :: symbol
    real(wp),intent(in),optional :: paramar(11)
    integer :: last_pg,i

    ! Reset state
    call reset_state()

    ! Set parameters if provided
    if (present(paramar)) then
      verbose = nint(paramar(1))
      MaxAxisOrder = nint(paramar(2))
      MaxOptCycles = nint(paramar(3))
      ToleranceSame = paramar(4)
      TolerancePrimary = paramar(5)
      ToleranceFinal = paramar(6)
      MaxOptStep = paramar(7)
      MinOptStep = paramar(8)
      GradientStep = paramar(9)
      OptChangeThreshold = paramar(10)
      OptChangeHits = nint(paramar(11))
    end if

    ! Set up atoms
    AtomsCount = natoms
    allocate (Atoms(AtomsCount))

    do i = 1,AtomsCount
      Atoms(i)%atom_type = attype(i)
      Atoms(i)%x(1) = coord(1,i)
      Atoms(i)%x(2) = coord(2,i)
      Atoms(i)%x(3) = coord(3,i)
    end do

    ! Find and analyze symmetry
    call find_symmetry_elements()
    call sort_symmetry_elements()
    call summarize_symmetry_elements()
    call report_symmetry_elements_brief()

    last_pg = identify_point_group()

    if (last_pg >= 1) then
      symbol = trim(PointGroups(last_pg)%group_name)
    else
      call report_symmetry_elements_brief_conly()
      if (len_trim(MaxRotAxis) == 0) then
        symbol = "C1"
      else
        symbol = trim(MaxRotAxis)
      end if
    end if
  end subroutine schoenflies

! ══════════════════════════════════════════════════════════════════════════════
! ══════════════════════════════════════════════════════════════════════════════
end module symmetry_i
