
module irmsd_module
!*****************************************
!* Module that implements a more
!* modern interface to calculating RMSDs
!*****************************************
  use crest_parameters
  use ls_rmsd,only:rmsd_classic => rmsd
  use strucrd
  use hungarian_module
  implicit none
  private

  public :: rmsd

  real(wp),parameter :: bigval = huge(bigval)

  public :: rmsd_cache
  type :: rmsd_cache
!****************************************************
!* cache implementation to avoid repeated allocation
!* and enable shared-memory parallelism
!****************************************************
    real(wp),allocatable :: xyzscratch(:,:,:)
    integer,allocatable  :: rankscratch(:,:)
    integer,allocatable  :: orderscratch(:)
    logical,allocatable  :: assignedscratch(:)

    type(assignment_cache),allocatable :: acache
  contains
    procedure :: allocate => allocate_rmsd_cache
  end type rmsd_cache

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine allocate_rmsd_cache(self,nat)
    implicit none
    class(rmsd_cache),intent(inout) :: self
    integer,intent(in) :: nat
    if (allocated(self%xyzscratch)) deallocate (self%xyzscratch)
    if (allocated(self%rankscratch)) deallocate (self%rankscratch)
    if (allocated(self%orderscratch)) deallocate (self%orderscratch)
    if (allocated(self%assignedscratch)) deallocate (self%assignedscratch)
    if (allocated(self%acache)) deallocate (self%acache)
    allocate (self%assignedscratch(nat),source=.false.)
    allocate (self%orderscratch(nat),source=0)
    allocate (self%rankscratch(nat,2),source=0)
    allocate (self%xyzscratch(3,nat,2),source=0.0_wp)
    allocate (self%acache)
    call self%acache%allocate(nat,nat)
  end subroutine allocate_rmsd_cache

  function rmsd(ref,mol,mask,scratch,rotmat,gradient) result(rmsdval)
!************************************************************************
!* function rmsd
!* Calculate the molecular RMSD via a quaternion algorithm
!*
!* Optional arguments are
!*   mask - boolean array to select a substructure for RMSD calculation
!*   scratch - workspace to create the substructures
!*   rotmat  - rotation matrix as return argument
!*   gradient - Cartesian gradient of the RMSD
!************************************************************************
    implicit none
    real(wp) :: rmsdval
    type(coord),intent(in) :: ref
    type(coord),intent(in) :: mol
    !> OPTIONAL arguments
    logical,intent(in),optional            :: mask(ref%nat)
    real(wp),intent(inout),target,optional :: scratch(3,ref%nat,2)
    real(wp),intent(out),optional          :: rotmat(3,3)
    real(wp),intent(out),target,optional   :: gradient(3,ref%nat)
    !> variables
    real(wp) :: x_center(3),y_center(3),Udum(3,3)
    real(wp),target :: gdum(3,3)
    integer  :: nat,getrotmat
    real(wp),allocatable,target :: tmpscratch(:,:,:)
    logical :: getgrad
    real(wp),pointer :: grdptr(:,:)
    real(wp),pointer :: scratchptr(:,:,:)
    integer :: ic,k

    !> initialize to large value
    rmsdval = bigval
    !> check structure consistency
    if (mol%nat .ne. ref%nat) return

    !> get rotation matrix?
    getrotmat = 0
    if (present(rotmat)) getrotmat = 1

    !> get gradient?
    if (present(gradient)) then
      getgrad = .true.
      gradient(:,:) = 0.0_wp
      grdptr => gradient
    else
      getgrad = .false.
      grdptr => gdum
    end if

!>--- substructure?
    if (present(mask)) then
      nat = count(mask(:))
      !> scratch workspace to use?
      if (present(scratch)) then
        scratchptr => scratch
      else
        allocate (tmpscratch(3,nat,2))
        scratchptr => tmpscratch
      end if

      !> do the mapping
      k = 0
      do ic = 1,ref%nat
        if (mask(ic)) then
          k = k+1
          scratchptr(1:3,k,1) = mol%xyz(1:3,ic)
          scratchptr(1:3,k,2) = ref%xyz(1:3,ic)
        end if
      end do

      !> calculate
      call rmsd_classic(nat,scratchptr(1:3,1:nat,1),scratchptr(1:3,1:nat,2), &
      &    getrotmat,Udum,x_center,y_center,rmsdval, &
      &    getgrad,grdptr)

      !> go backwards through gradient (if necessary) to restore atom order
      if (getgrad) then
        k = nat
        do ic = nat,1,-1
          if (mask(ic)) then
            grdptr(1:3,ic) = grdptr(1:3,k)
            grdptr(1:3,k) = 0.0_wp
            k = k-1
          end if
        end do
      end if

      nullify(scratchptr)
      if (allocated(tmpscratch)) deallocate (tmpscratch)

    else
!>--- standard calculation (Quarternion algorithm)
      call rmsd_classic(ref%nat,mol%xyz,ref%xyz, &
      &    getrotmat,Udum,x_center,y_center,rmsdval, &
      &    getgrad,grdptr)
    end if

    !> pass on rotation matrix if asked for
    if (getrotmat > 0) rotmat = Udum

  end function rmsd

!========================================================================================!

  subroutine min_rmsd(ref,mol,rcache,rmsdout)
    implicit none
    !> IN & OUTPUT
    type(coord),intent(in)    :: ref
    type(coord),intent(inout) :: mol
    type(rmsd_cache),intent(inout),optional,target :: rcache
    real(wp),intent(out),optional :: rmsdout

    !> LOCAL
    type(rmsd_cache),pointer :: cptr
    type(rmsd_cache),allocatable,target :: local_rcache
    integer :: natmax
    real(wp) :: calc_rmsd

    if (present(rcache)) then
      cptr => rcache
    else
      allocate (local_rcache)
      natmax = max(ref%nat,mol%nat)
      call local_rcache%allocate(natmax)
      cptr => local_rcache
    end if

    calc_rmsd = rmsd(ref,mol,scratch=cptr%xyzscratch)

    if (present(rmsdout)) rmsdout = calc_rmsd
  end subroutine min_rmsd

!========================================================================================!

  subroutine atswp(mol,ati,atj)
  !********************************
  !* swap atom ati with atj in mol
  !********************************
     implicit none
     type(coord),intent(inout) :: mol
     integer,intent(in) :: ati,atj 
     real(wp) :: xyztmp(3)
     integer :: attmp
     xyztmp(1:3) = mol%xyz(1:3,ati)
     attmp = mol%at(ati)
     mol%xyz(1:3,ati) = mol%xyz(1:3,atj)
     mol%at(ati) = mol%at(atj)
     mol%xyz(1:3,atj) = xyztmp(1:3)
     mol%at(atj) = attmp
  end subroutine atswp

!========================================================================================!

  subroutine compute_hungarian(ref,mol,acache)
    implicit none
    !> IN & OUTPUT
    type(coord),intent(in)    :: ref 
    type(coord),intent(inout) :: mol
    type(assignment_cache),intent(inout),optional,target :: acache

    !> LOCAL
    type(assignment_cache),pointer :: aptr
    type(assignment_cache),allocatable,target :: local_acache
    integer :: natmax
    

    if (present(acache)) then
      aptr => acache
    else
      allocate (local_acache)
      natmax = max(ref%nat,mol%nat)
      call local_acache%allocate(natmax,natmax)
      aptr => local_acache
    end if

    !> Compute the cost matrix, which is simply the distance matrix 
    !> between the two molecules.
    !> To avoid computational overhead we can skip the square root. 
    !> It won't affect the result 



  end subroutine compute_hungarian

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
end module irmsd_module
