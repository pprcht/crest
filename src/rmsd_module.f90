

module rmsd_module
!*****************************************
!* Module that implements a more
!* modern interface to calculating RMSDs
!*****************************************
    use crest_parameters
    use ls_rmsd, only: rmsd_classic => rmsd
    use strucrd
    implicit none
    private 

    public :: rmsd


    real(wp),parameter :: bigval = huge(bigval)

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

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
    real(wp) :: x_center(3),y_center(3),Udum(3,3),gdum(3,3)
    integer  :: nat,getrotmat 
    real(wp),allocatable :: tmpscratch(:,:,:)
    logical :: getgrad
    real(wp),pointer :: grdptr !(:,:)
    real(wp),pointer :: scratchptr !(:,:,:)
    integer :: ic,k

    !> initialize to large value
    rmsdval = bigval
    !> check structure consistency
    if(mol%nat .ne. ref%nat) return

    !> get rotation matrix?
    getrotmat = 0
    if(present(rotmat)) getrotmat = 1

    !> get gradient?
    if(present(gradient))then
      getgrad = .true.
      gradient(:,:) = 0.0_wp
      grdptr => gradient
    else
      getgrad = .false.
      grdptr => gdum
    endif 

!>--- substructure?
    if(present(mask))then
      nat = count(mask(:))
      !> scratch workspace to use?
      if(present(scratch))then
        scratchptr => scratch
      else
        allocate(tmpscratch(3,nat,2))
        scratchptr => tmpscratch
      endif

      !> do the mapping
      k=0
      do ic=1,ref%nat
        if(mask(ic))then
          k=k+1
          scratchptr(1:3,k,1) = mol%xyz(1:3,ic)
          scratchptr(1:3,k,2) = ref%xyz(1:3,ic) 
        endif
      enddo

      !> calculate
      call rmsd_classic(nat, scratchptr(1:3,1:nat,1), scratchptr(1:3,1:nat,2), &
      &    getrotmat, Udum, x_center, y_center, rmsdval, &
      &    getgrad, grdptr)

      !> go backwards through gradient (if necessary) to restore atom order
      if(getgrad)then
        k=nat
        do ic=nat,1,-1 
          if(mask(ic))then
            grdptr(1:3,ic) = grdptr(1:3,k)
            grdptr(1:3,k) = 0.0_wp
            k=k-1
          endif 
        enddo
      endif
    
      deallocate(scratchptr)
      if(allocated(tmpscratch)) deallocate(tmpscratch)

    else
!>--- standard calculation
      call rmsd_classic(ref%nat,mol%xyz,ref%xyz, &
      &    getrotmat, Udum, x_center, y_center, rmsdval, &
      &    getgrad, grdptr)
    endif

    !> pass on rotation matrix if asked for
    if(getrotmat > 0) rotmat = Udum

end function rmsd

!========================================================================================!
!========================================================================================!
end module rmsd_module
