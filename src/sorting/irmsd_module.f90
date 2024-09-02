
module irmsd_module
!*****************************************
!* Module that implements a more
!* modern interface to calculating RMSDs
!*****************************************
  use crest_parameters
  use ls_rmsd,only:rmsd_classic => rmsd
  use strucrd
  use hungarian_module
  use axis_module
  implicit none
  private

  public :: rmsd
  public :: min_rmsd

  real(wp),parameter :: bigval = huge(bigval)

  public :: rmsd_cache
  type :: rmsd_cache
!****************************************************
!* cache implementation to avoid repeated allocation
!* and enable shared-memory parallelism
!****************************************************
    real(wp),allocatable :: xyzscratch(:,:,:)
    integer,allocatable  :: rank(:,:)
    integer,allocatable  :: best_order(:,:)
    integer,allocatable  :: current_order(:)
    integer,allocatable  :: target_order(:)
    integer,allocatable  :: iwork(:)
    logical,allocatable  :: assigned(:)  !> atom-wise
    logical,allocatable  :: rassigned(:) !> rank-wise

    integer :: nranks = 0
    integer,allocatable :: ngroup(:)
    logical :: stereocheck = .false.

    type(assignment_cache),allocatable :: acache
  contains
    procedure :: allocate => allocate_rmsd_cache
  end type rmsd_cache

  real(wp),parameter :: inf = huge(1.0_wp)

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
    if (allocated(self%rank)) deallocate (self%rank)
    if (allocated(self%best_order)) deallocate (self%best_order)
    if (allocated(self%current_order)) deallocate (self%current_order)
    if (allocated(self%target_order)) deallocate (self%target_order)
    if (allocated(self%iwork)) deallocate (self%iwork)
    if (allocated(self%assigned)) deallocate (self%assigned)
    if (allocated(self%rassigned)) deallocate (self%rassigned)
    if (allocated(self%ngroup)) deallocate (self%ngroup)
    if (allocated(self%acache)) deallocate (self%acache)
    allocate (self%assigned(nat),source=.false.)
    allocate (self%rassigned(nat),source=.false.)
    allocate (self%best_order(nat,3),source=0)
    allocate (self%current_order(nat),source=0)
    allocate (self%target_order(nat),source=0)
    allocate (self%iwork(nat),source=0)
    allocate (self%rank(nat,2),source=0)
    self%nranks = 0
    allocate (self%ngroup(nat), source=0)
    allocate (self%xyzscratch(3,nat,2),source=0.0_wp)
    allocate (self%acache)
    call self%acache%allocate(nat,nat,.true.) !> assume we are only using the LSAP implementation
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

      nullify (scratchptr)
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
    integer :: nat,ii,rnk
    real(wp) :: calc_rmsd
    real(wp) :: tmprmsd_sym(3)
    logical,parameter :: debug = .true.

!>--- Initialization
    if (present(rcache)) then
      cptr => rcache
    else
      allocate (local_rcache)
      if (ref%nat .ne. mol%nat) then
        error stop 'Unequal molecule size in min_rmsd()'
      end if
      nat = max(ref%nat,mol%nat)
      call local_rcache%allocate(nat)
      call fallbackranks(ref,mol,nat,local_rcache%rank)
      cptr => local_rcache
    end if

!>-- Consistency check
    cptr%nranks = maxval(cptr%rank(:,1))
    if(cptr%nranks .ne. maxval(cptr%rank(:,2)))then
      error stop "Different atom identities in min_rmsd, can't restore an atom order!"
    endif

!>--- First sorting, to at least restore rank order
    if (.not.all(cptr%rank(:,1) .eq. cptr%rank(:,2))) then
      call rank_2_order(ref%nat,cptr%rank(:,1),cptr%target_order)
      call rank_2_order(mol%nat,cptr%rank(:,2),cptr%current_order)
      if (debug) then
        write (*,*) 'current order & rank & target order'
        do ii = 1,mol%nat
          write (*,*) cptr%current_order(ii),cptr%rank(ii,2),cptr%target_order(ii)
        end do
      end if
      call molatomsort(mol,mol%nat,cptr%current_order,cptr%target_order,cptr%iwork)
      cptr%rank(:,2) = cptr%rank(:,1) !> since the ranks must be equal now!
      if (debug) then
        write (*,*) 'sorted order & rank'
        do ii = 1,mol%nat
          write (*,*) cptr%current_order(ii),cptr%rank(ii,2)
        end do
      end if
    end if

!>--- Count symmetry equivalent groups and assign all unique atoms immediately   
!     Note, the rank can be zero if we only are looking at heavy atoms
    if(all(cptr%ngroup(:) .eq. 0))then 
      do ii=1,ref%nat
         rnk = cptr%rank(ii,1)
         if(rnk>0)then
           cptr%ngroup(rnk) = cptr%ngroup(rnk) + 1
         endif
      enddo
    endif
    !> assignment reset
    cptr%assigned(:) = .false.
    cptr%rassigned(:) = .false.
    cptr%rassigned(cptr%nranks:) = .true.
    do ii=1,ref%nat
      rnk = cptr%rank(ii,2)
      if(rnk < 1)then
        cptr%assigned(ii) = .true.
        cycle
      endif
      if(cptr%ngroup(rnk) .eq. 1)then
        cptr%assigned(ii) = .true.
      endif
    enddo
    if(debug)then
      write (*,*)  'rank & # members'
      do ii = 1,mol%nat
        if(cptr%ngroup(ii) > 0)then
          write (*,*) ii,cptr%ngroup(ii)
        endif
      end do
    endif

!>--- Perform the desired symmetry operations, align with rotational axis, run LSAP algo
    tmprmsd_sym(:) = inf  !> initialize to huge
    call axis(mol%nat, mol%at, mol%xyz)
    call min_rmsd_iterate_through_groups(ref,mol,1,cptr)

    !> mirror z
    if(cptr%stereocheck)then
      mol%xyz(3,:) = -mol%xyz(3,:)
      call axis(mol%nat, mol%at, mol%xyz)
      call min_rmsd_iterate_through_groups(ref,mol,1,cptr)
    endif


!>--- select the best match among the ones after symmetry operations and use its ordering


!>--- final RMSD with fully restored atom order
    calc_rmsd = rmsd(ref,mol,scratch=cptr%xyzscratch)

    if (present(rmsdout)) rmsdout = calc_rmsd
  end subroutine min_rmsd

!========================================================================================!

  subroutine min_rmsd_iterate_through_groups(ref,mol,step,rcache)
      implicit none
      type(coord),intent(in) :: ref
      type(coord),intent(inout) :: mol
      integer,intent(in) :: step
      type(rmsd_cache),intent(inout) :: rcache
      integer :: rr
      logical,parameter :: debug=.true.

      if(debug)then
      write(*,*) '# ranks:',rcache%nranks
      endif
      do rr=1,rcache%nranks
        if(rcache%rassigned(rr)) cycle
 
      enddo

  end subroutine min_rmsd_iterate_through_groups

!========================================================================================!

  subroutine fallbackranks(ref,mol,nat,ranks)
!*****************************************************************
!* If we are doing ranks on-the-fly (i.e. without canonical algo)
!* we can fall back to just using the atom types
!*****************************************************************
    implicit none
    type(coord),intent(in) :: ref,mol
    integer,intent(in)     :: nat
    integer,intent(inout)  :: ranks(nat,2)

    integer,allocatable :: typemap(:),rtypemap(:)
    integer :: k,ii
    allocate (typemap(nat),source=0)
    k = 0
    do ii = 1,ref%nat
      if (.not.any(typemap(:) .eq. ref%at(ii))) then
        k = k+1
        typemap(k) = ref%at(ii)
      end if
    end do
    do ii = 1,mol%nat
      if (.not.any(typemap(:) .eq. mol%at(ii))) then
        k = k+1
        typemap(k) = mol%at(ii)
      end if
    end do
    k = maxval(typemap(:))
    allocate (rtypemap(k),source=0)
    do ii = 1,nat
      if (typemap(ii) == 0) cycle
      rtypemap(typemap(ii)) = ii
    end do
    !> assign
    do ii = 1,ref%nat
      ranks(ii,1) = rtypemap(ref%at(ii))
    end do
    do ii = 1,mol%nat
      ranks(ii,2) = rtypemap(mol%at(ii))
    end do
    deallocate (rtypemap)
    deallocate (typemap)
  end subroutine fallbackranks

!========================================================================================!

  subroutine compute_hungarian(ref,mol,ranks,targetrank,acache)
!**************************************************************
!* Run the linear assignment algorithm on the desired subset
!* of atoms (via rank and targetrank)
!**************************************************************
    implicit none
    !> IN & OUTPUT
    type(coord),intent(in)    :: ref
    type(coord),intent(inout) :: mol
    integer,intent(in) :: ranks(:)
    integer,intent(in) :: targetrank
    type(assignment_cache),intent(inout),optional,target :: acache

    !> LOCAL
    type(assignment_cache),pointer :: aptr
    type(assignment_cache),allocatable,target :: local_acache
    integer :: nat,i,ii,jj

    if (present(acache)) then
      aptr => acache
    else
      allocate (local_acache)
      if (ref%nat .ne. mol%nat) then
        error stop 'Unequal molecule size in compute_hungarian()'
      end if
      nat = max(ref%nat,mol%nat)
      call local_acache%allocate(nat,nat,.true.)
      aptr => local_acache
    end if

    !> Compute the cost matrix, which is simply the distance matrix
    !> between the two molecules.
    !> To avoid computational overhead we can skip the square root.
    !> It won't affect the result

  end subroutine compute_hungarian


!========================================================================================!



!========================================================================================!

  subroutine rank_2_order(nat,rank,order)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: rank(nat)
    integer,intent(out) :: order(nat)
    integer :: ii,jj,k,maxrank
    order(:) = 0
    maxrank = maxval(rank(:))
    k = 0
    do ii = 1,maxrank
      do jj = 1,nat
        if (rank(jj) == ii) then
          k = k+1
          order(jj) = k
        end if
      end do
    end do
  end subroutine rank_2_order

!========================================================================================!

  subroutine molatomsort(mol,n,current_order,target_order,index_map)
    implicit none
    type(coord),intent(inout) :: mol
    integer,intent(in) :: n
    integer,intent(inout) :: current_order(n)
    integer,intent(in) :: target_order(n)
    integer,intent(inout) :: index_map(n)
    integer :: i,j,correct_atom,current_position

    !> Step 1: Create a mapping from target_order to current_order positions
    do i = 1,n
      index_map(current_order(i)) = i
    end do

    !> Step 2: Restore the target order
    do i = 1,n
      correct_atom = target_order(i)
      current_position = index_map(correct_atom)

      if (i /= current_position) then
        !> Swap atoms i and current_position in molecule
        call mol%swap(i,current_position)

        !> Update the index map since the atoms have been swapped
        index_map(current_order(i)) = current_position
        index_map(current_order(current_position)) = i

        !> Update the current_order array to reflect the swap
        j = current_order(i)
        current_order(i) = current_order(current_position)
        current_order(current_position) = j
      end if
    end do
  end subroutine molatomsort

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
end module irmsd_module
