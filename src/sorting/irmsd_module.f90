
module irmsd_module
!*****************************************
!* Module that implements a more
!* modern interface to calculating RMSDs
!*****************************************
  use crest_parameters
  use strucrd
  use hungarian_module
  use axis_module
  implicit none
  private

  public :: rmsd
  public :: min_rmsd

  public :: checkranks,fallbackranks

  real(wp),parameter :: bigval = huge(bigval)

  type :: rmsd_core_cache
!*************************************
!* Memory cache for rmsd_core routine
!*************************************
    real(wp),allocatable :: x(:,:)
    real(wp),allocatable :: y(:,:)
    real(wp),allocatable :: xi(:)
    real(wp),allocatable :: yi(:)
  contains
    procedure :: allocate => allocate_rmsd_core_cache
  end type rmsd_core_cache

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
    integer,allocatable  :: order_bkup(:,:)
    integer,allocatable  :: iwork(:)
    integer,allocatable  :: iwork2(:,:)
    logical,allocatable  :: assigned(:)  !> atom-wise
    logical,allocatable  :: rassigned(:) !> rank-wise

    integer :: nranks = 0
    integer,allocatable :: ngroup(:)
    logical :: stereocheck = .false.

    type(rmsd_core_cache),allocatable :: ccache
    type(assignment_cache),allocatable :: acache
  contains
    procedure :: allocate => allocate_rmsd_cache
  end type rmsd_cache

  real(wp),parameter :: inf = huge(1.0_wp)
  real(wp),parameter :: imat(3,3) = reshape([1.0_wp,0.0_wp,0.0_wp, &
                                    &        0.0_wp,1.0_wp,0.0_wp, &
                                    &        0.0_wp,0.0_wp,1.0_wp], &
                                    &        [3,3])

  real(wp),parameter :: Rx180(3,3) = reshape([1.0_wp,0.0_wp,0.0_wp,   &
                                     &          0.0_wp,-1.0_wp,0.0_wp,   &
                                     &          0.0_wp,0.0_wp,-1.0_wp], &
                                     &          [3,3])

  real(wp),parameter :: Ry180(3,3) = reshape([-1.0_wp,0.0_wp,0.0_wp,   &
                                     &          0.0_wp,1.0_wp,0.0_wp,   &
                                     &          0.0_wp,0.0_wp,-1.0_wp], &
                                     &          [3,3])

  real(wp),parameter :: Rz180(3,3) = reshape([-1.0_wp,0.0_wp,0.0_wp,   &
                                     &          0.0_wp,-1.0_wp,0.0_wp,   &
                                     &          0.0_wp,0.0_wp,1.0_wp], &
                                     &          [3,3])

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine allocate_rmsd_core_cache(self,nat)
    implicit none
    class(rmsd_core_cache),intent(inout) :: self
    integer,intent(in) :: nat
    if (allocated(self%x)) deallocate (self%x)
    if (allocated(self%y)) deallocate (self%y)
    if (allocated(self%xi)) deallocate (self%xi)
    if (allocated(self%yi)) deallocate (self%yi)
    allocate (self%xi(nat),source=0.0_wp)
    allocate (self%yi(nat),source=0.0_wp)
    allocate (self%x(3,nat),source=0.0_wp)
    allocate (self%y(3,nat),source=0.0_wp)
  end subroutine allocate_rmsd_core_cache

  subroutine allocate_rmsd_cache(self,nat)
    implicit none
    class(rmsd_cache),intent(inout) :: self
    integer,intent(in) :: nat
    if (allocated(self%xyzscratch)) deallocate (self%xyzscratch)
    if (allocated(self%rank)) deallocate (self%rank)
    if (allocated(self%best_order)) deallocate (self%best_order)
    if (allocated(self%current_order)) deallocate (self%current_order)
    if (allocated(self%target_order)) deallocate (self%target_order)
    if (allocated(self%order_bkup)) deallocate (self%order_bkup)
    if (allocated(self%iwork)) deallocate (self%iwork)
    if (allocated(self%iwork2)) deallocate (self%iwork2)
    if (allocated(self%assigned)) deallocate (self%assigned)
    if (allocated(self%rassigned)) deallocate (self%rassigned)
    if (allocated(self%ngroup)) deallocate (self%ngroup)
    if (allocated(self%ccache)) deallocate (self%ccache)
    if (allocated(self%acache)) deallocate (self%acache)
    allocate (self%assigned(nat),source=.false.)
    allocate (self%rassigned(nat),source=.false.)
    allocate (self%best_order(nat,3),source=0)
    allocate (self%current_order(nat),source=0)
    allocate (self%target_order(nat),source=0)
    allocate (self%order_bkup(nat,8),source=0)
    allocate (self%iwork(nat),source=0)
    allocate (self%iwork2(nat,2),source=0)
    allocate (self%rank(nat,2),source=0)
    self%nranks = 0
    allocate (self%ngroup(nat),source=0)
    allocate (self%xyzscratch(3,nat,2),source=0.0_wp)
    allocate (self%ccache)
    allocate (self%acache)
    call self%ccache%allocate(nat)
    call self%acache%allocate(nat,nat,.true.) !> assume we are only using the LSAP implementation
  end subroutine allocate_rmsd_cache

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

  function rmsd(ref,mol,mask,scratch,rotmat,gradient,ccache) result(rmsdval)
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
    type(rmsd_core_cache),intent(inout),optional,target :: ccache
    !> variables
    type(rmsd_core_cache),allocatable,target :: ccachetmp
    type(rmsd_core_cache),pointer :: ccptr
    real(wp) :: x_center(3),y_center(3),Udum(3,3)
    real(wp),target :: gdum(3,3)
    integer  :: nat,getrotmat
    logical  :: calc_u
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
    calc_u = .false.
    if (present(rotmat)) then
      getrotmat = 1
      calc_u = .true.
    end if

    !> get gradient?
    if (present(gradient)) then
      getgrad = .true.
      gradient(:,:) = 0.0_wp
      grdptr => gradient
    else
      getgrad = .false.
      grdptr => gdum
    end if

    !> use present cache?
    if (present(ccache)) then
      ccptr => ccache
    else
      allocate (ccachetmp)
      call ccachetmp%allocate(ref%nat)
      ccptr => ccachetmp
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
      call rmsd_core(nat,scratchptr(1:3,1:nat,1),scratchptr(1:3,1:nat,2), &
      &          calc_u,Udum,rmsdval,getgrad,grdptr,ccptr)

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
      call rmsd_core(ref%nat,mol%xyz,ref%xyz, &
      &          calc_u,Udum,rmsdval,getgrad,grdptr,ccptr)
    end if

    !> pass on rotation matrix if asked for
    if (calc_u) rotmat = Udum

  end function rmsd

!========================================================================================!

  subroutine rmsd_core(nat,xyz1,xyz2,calc_u,U,error,calc_g,grad,ccache)
!**********************************************************
!* Rewrite or RMSD code with modified memory management
!* Adapted from ls_rmsd, and using some of its subroutines
!* The goal is to offload memory allocation to outside
!* the routine in case it is repeadetly called
!**********************************************************
    use ls_rmsd,only:dstmev,rotation_matrix
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: xyz1(3,nat)
    real(wp),intent(in) :: xyz2(3,nat)
    logical,intent(in) :: calc_u
    real(wp),dimension(3,3),intent(out) :: U
    real(wp),intent(out) :: error
    logical,intent(in) :: calc_g
    real(wp),intent(inout) :: grad(:,:)
    type(rmsd_core_cache),intent(inout) :: ccache

    !> LOCAL
    integer :: i,j
    real(wp) :: x_center(3)
    real(wp) :: y_center(3)
    real(wp) :: x_norm,y_norm,lambda
    real(wp) :: Rmatrix(3,3)
    real(wp) :: S(4,4)
    real(wp) :: q(4)
    real(wp) :: tmp(3),rnat
    integer :: io

    !> associate
    associate (x => ccache%x,y => ccache%y,xi => ccache%xi,yi => ccache%yi)

      !> make copies of the original coordinates
      x(:,:) = xyz1(:,:)
      y(:,:) = xyz2(:,:)

      !> calculate the barycenters, centroidal coordinates, and the norms
      x_norm = 0.0_wp
      y_norm = 0.0_wp
      rnat = 1.0_wp/real(nat,wp)
      do i = 1,3
        xi(:) = x(i,:)
        yi(:) = y(i,:)
        x_center(i) = sum(xi(1:nat))*rnat
        y_center(i) = sum(yi(1:nat))*rnat
        xi(:) = xi(:)-x_center(i)
        yi(:) = yi(:)-y_center(i)
        x(i,:) = xi(:)
        y(i,:) = yi(:)
        x_norm = x_norm+dot_product(xi,xi)
        y_norm = y_norm+dot_product(yi,yi)
      end do

      !> calculate the R matrix
      do i = 1,3
        do j = 1,3
          Rmatrix(i,j) = dot_product(x(i,:),y(j,:))
        end do
      end do

      !> S matrix
      S(1,1) = Rmatrix(1,1)+Rmatrix(2,2)+Rmatrix(3,3)
      S(2,1) = Rmatrix(2,3)-Rmatrix(3,2)
      S(3,1) = Rmatrix(3,1)-Rmatrix(1,3)
      S(4,1) = Rmatrix(1,2)-Rmatrix(2,1)

      S(1,2) = S(2,1)
      S(2,2) = Rmatrix(1,1)-Rmatrix(2,2)-Rmatrix(3,3)
      S(3,2) = Rmatrix(1,2)+Rmatrix(2,1)
      S(4,2) = Rmatrix(1,3)+Rmatrix(3,1)

      S(1,3) = S(3,1)
      S(2,3) = S(3,2)
      S(3,3) = -Rmatrix(1,1)+Rmatrix(2,2)-Rmatrix(3,3)
      S(4,3) = Rmatrix(2,3)+Rmatrix(3,2)

      S(1,4) = S(4,1)
      S(2,4) = S(4,2)
      S(3,4) = S(4,3)
      S(4,4) = -Rmatrix(1,1)-Rmatrix(2,2)+Rmatrix(3,3)

      !> Calculate eigenvalues and eigenvectors, and
      !> take the maximum eigenvalue lambda and the corresponding eigenvector q.
      call dstmev(S,lambda,q,io)
      if (io /= 0) then
        error = -1.0_wp
        return
      end if

      if (calc_u) then
        !> reset
        U(:,:) = Imat(:,:)
        !> convert quaternion q to rotation matrix U
        call rotation_matrix(q,U)
      end if

      !> RMS Deviation
      error = sqrt(max(0.0_wp, ((x_norm+y_norm)-2.0_wp*lambda))*rnat)

      if (calc_g) then
        !> Gradient of the error of xyz1 w.r.t xyz2
        do i = 1,nat
          do j = 1,3
            tmp(:) = matmul(transpose(U(:,:)),y(:,i))
            grad(j,i) = ((x(j,i)-tmp(j))/error)*rnat
          end do
        end do
      end if

    end associate
  end subroutine rmsd_core

!========================================================================================!

  subroutine min_rmsd(ref,mol,rcache,rmsdout,align)
!*********************************************************************
!* Main routine to determine minium RMSD considering atom permutation
!* Input
!*   ref  - the reference structure
!*   mol  - the structure to be matched to ref
!* Optinal arguments
!*   rcache  - memory cache
!*   rmsdout - the calculated RMSD scalar
!*   align   - quarternion-align mol in the last stage
!*********************************************************************
    implicit none
    !> IN & OUTPUT
    type(coord),intent(in)    :: ref
    type(coord),intent(inout) :: mol
    type(rmsd_cache),intent(inout),optional,target :: rcache
    real(wp),intent(out),optional :: rmsdout
    logical,intent(in),optional :: align

    !> LOCAL
    type(rmsd_cache),pointer :: cptr
    type(rmsd_cache),allocatable,target :: local_rcache
    integer :: nat,ii,rnk,dumpunit
    real(wp) :: calc_rmsd
    real(wp) :: tmprmsd_sym(8),dum
    real(wp) :: rotmat(3,3)
    logical,parameter :: debug = .false.

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
    if (cptr%nranks .ne. maxval(cptr%rank(:,2))) then
      error stop "Different atom identities in min_rmsd, can't restore an atom order!"
    end if

!>--- First sorting, to at least restore rank order (only if that's not the case!)
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
    if (all(cptr%ngroup(:) .eq. 0)) then
      do ii = 1,ref%nat
        rnk = cptr%rank(ii,1)
        if (rnk > 0) then
          cptr%ngroup(rnk) = cptr%ngroup(rnk)+1
        end if
      end do
    end if
    !> assignment reset
    cptr%assigned(:) = .false.
    cptr%rassigned(:) = .false.
    cptr%rassigned(cptr%nranks+1:) = .true. !> skip unneeded allocation space
    do ii = 1,ref%nat
      cptr%iwork(ii) = ii         !> also init iwork
      cptr%target_order(ii) = ii  !> also init target_order
      rnk = cptr%rank(ii,2)
      if (rnk < 1) then
        cptr%assigned(ii) = .true.
        cycle
      end if
      if (cptr%ngroup(rnk) .eq. 1) then
        cptr%assigned(ii) = .true.
        cptr%rassigned(rnk) = .true.
      end if
    end do
    if (debug) then
      write (*,*) 'rank & # members'
      do ii = 1,mol%nat
        if (cptr%ngroup(ii) > 0) then
          write (*,*) ii,cptr%ngroup(ii)
        end if
      end do
    end if

!>--- Perform the desired symmetry operations, align with rotational axis, run LSAP algo
!>    Since the rotational axis alignment can be a bit arbitrary w.r.t 180° rotations
!>    we need to check these as well.
    if (debug) then
      open (newunit=dumpunit,file='debugirmsd.xyz')
      call ref%append(dumpunit)
    end if
    !> initialize to huge
    tmprmsd_sym(:) = inf
    !> initial alignment of mol
    call axis(mol%nat,mol%at,mol%xyz)

    !> Running the checks
    call min_rmsd_rotcheck(ref,mol,cptr,tmprmsd_sym,1)
    if (debug) then
      write (*,*) 'Total LSAP cost:',minval(tmprmsd_sym(1:4))
      call mol%append(dumpunit)
    end if

    !> mirror z and re-run the same checks
    if (cptr%stereocheck) then
      mol%xyz(3,:) = -mol%xyz(3,:)  !> mirror z
      call axis(mol%nat,mol%at,mol%xyz) !> align

      !> Running the checks
      call min_rmsd_rotcheck(ref,mol,cptr,tmprmsd_sym,2)
      if (debug) then
        write (*,*) 'Total LSAP cost (inverted):',minval(tmprmsd_sym(5:8))
        call mol%append(dumpunit)
      end if
      mol%xyz(3,:) = -mol%xyz(3,:)  !> restore z
    end if

!>--- select the best match among the ones after symmetry operations and use its ordering
    ii = minloc(tmprmsd_sym(:),1)
    if (debug) then
      write (*,*) 'final alignment:',ii,"/ 8"
    end if
    if (ii > 4) then
      mol%xyz(3,:) = -mol%xyz(3,:)
      if (debug) write (*,*) 'inverting'
    end if
    select case (ii) !> 180° rotations
    case (1,5)
      continue
    case (2,6)
      mol%xyz = matmul(Rx180,mol%xyz)
      if (debug) write (*,*) '180°x'
    case (3,7)
      mol%xyz = matmul(Rx180,mol%xyz)
      mol%xyz = matmul(Ry180,mol%xyz)
      if (debug) write (*,*) '180°x, 180°y'
    case (4,8)
      mol%xyz = matmul(Ry180,mol%xyz)
      if (debug) write (*,*) '180°y'
    end select
    cptr%current_order(:) = cptr%order_bkup(:,ii)

    if (debug) then
      write (*,*) 'Determined remapping'
      do ii = 1,mol%nat
        write (*,*) cptr%current_order(ii),'-->',cptr%target_order(ii)
      end do
    end if

    call molatomsort(mol,mol%nat,cptr%current_order,cptr%target_order,cptr%iwork)
    if (debug) then
      call mol%append(dumpunit)
      close (dumpunit)
    end if

!>--- final RMSD with fully restored atom order
    if (present(align)) then
      calc_rmsd = rmsd(ref,mol,scratch=cptr%xyzscratch,ccache=cptr%ccache,rotmat=rotmat)
      if (align) then
        mol%xyz = matmul(rotmat,mol%xyz)
      end if
    else
      calc_rmsd = rmsd(ref,mol,scratch=cptr%xyzscratch,ccache=cptr%ccache)
    end if

    if (present(rmsdout)) rmsdout = calc_rmsd
  end subroutine min_rmsd

!========================================================================================!

  subroutine min_rmsd_iterate_through_groups(ref,mol,rcache,val)
    implicit none
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    type(rmsd_cache),intent(inout),target :: rcache
    real(wp),intent(out) :: val
    integer :: rr,ii,jj
    real(wp) :: val0
    type(assignment_cache),pointer :: aptr
    logical,parameter :: debug = .false.

    !> reset val
    val = 0.0_wp

    if (debug) then
      write (*,*) '# ranks:',rcache%nranks
    end if
    aptr => rcache%acache
    do rr = 1,rcache%nranks
      if (rcache%rassigned(rr)) cycle

      !> LSAP wrapper that computes the relevant Cost matrix for the atoms of rank rr
      call compute_linear_sum_assignment( &
      &        ref,mol,rcache%rank,rcache%ngroup,rr, &
      &        rcache%iwork2,aptr,val0)

      do ii = 1,rcache%ngroup(rr)
        rcache%iwork(rcache%iwork2(ii,1)) = rcache%iwork2(ii,2)
      end do

      !> add up the total LSAP cost (of considered ranks)
      !> we need this if we have to decide on a mapping in case of false enantiomers
      val = val+val0
    end do

  end subroutine min_rmsd_iterate_through_groups

!========================================================================================!

  subroutine min_rmsd_rotcheck(ref,mol,cptr,values,step)
    implicit none
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    type(rmsd_cache),intent(inout),target :: cptr
    real(wp),intent(inout) :: values(:)
    integer,intent(in) :: step
    integer :: rr,ii,jj,debugunit2
    real(wp) :: vals(4),dum
    logical,parameter :: debug = .false.

    !> reset val
    vals(:) = 0.0_wp

    if (debug) then
      open (newunit=debugunit2,file='rotdebug.xyz')
      call ref%append(debugunit2)
    end if
    call min_rmsd_iterate_through_groups(ref,mol,cptr,dum)
    vals(1) = dum
    if (debug) call mol%append(debugunit2)
    cptr%order_bkup(:,1+4*(step-1)) = cptr%iwork(:)

    mol%xyz = matmul(Rx180,mol%xyz)
    call min_rmsd_iterate_through_groups(ref,mol,cptr,dum)
    vals(2) = dum
    if (debug) call mol%append(debugunit2)
    cptr%order_bkup(:,2+4*(step-1)) = cptr%iwork(:)

    mol%xyz = matmul(Ry180,mol%xyz)
    call min_rmsd_iterate_through_groups(ref,mol,cptr,dum)
    vals(3) = dum
    if (debug) call mol%append(debugunit2)
    cptr%order_bkup(:,3+4*(step-1)) = cptr%iwork(:)

    mol%xyz = matmul(Rx180,mol%xyz)
    call min_rmsd_iterate_through_groups(ref,mol,cptr,dum)
    vals(4) = dum
    if (debug) call mol%append(debugunit2)
    cptr%order_bkup(:,4+4*(step-1)) = cptr%iwork(:)

    mol%xyz = matmul(Ry180,mol%xyz) !> restore

    if (debug) then
      close (debugunit2)
      write (*,*) 'vals:',vals(:)
    end if

    do ii = 1,4
      values(ii+4*(step-1)) = vals(ii)
    end do
  end subroutine min_rmsd_rotcheck

!=========================================================================================!

  subroutine min_rmsd_quadalign(ref,mol,rotate)
    implicit none
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    logical,intent(out) :: rotate(3)
    integer :: rr,ii,jj,acount
    real(wp) :: val0
    integer :: tmp1,tmp2
    type(assignment_cache),pointer :: aptr
    logical,parameter :: debug = .false.

    rotate(:) = .false.

    do jj = 1,3
      tmp1 = count(ref%xyz(jj,:) > 0.0_wp)
      tmp2 = count(mol%xyz(jj,:) > 0.0_wp)
      if (tmp1 .ne. tmp2) rotate(jj) = .true.
    end do
  end subroutine min_rmsd_quadalign

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

  subroutine compute_linear_sum_assignment(ref,mol,ranks, &
                        & ngroups,targetrank,iwork2,acache,val0)
!**************************************************************
!* Run the linear assignment algorithm on the desired subset
!* of atoms (via rank and targetrank)
!**************************************************************
    implicit none
    !> IN & OUTPUT
    type(coord),intent(in)    :: ref
    type(coord),intent(inout) :: mol
    integer,intent(in) :: ranks(:,:)
    integer,intent(in) :: ngroups(:)
    integer,intent(in) :: targetrank
    integer,intent(inout) :: iwork2(:,:)
    type(assignment_cache),intent(inout),optional,target :: acache
    real(wp),intent(out) :: val0

    !> LOCAL
    type(assignment_cache),pointer :: aptr
    type(assignment_cache),allocatable,target :: local_acache
    integer :: nat,i,j,ii,jj,rnknat,iostatus
    real(sp) :: dists(3)

    logical,parameter :: debug = .false.

    val0 = 0.0_wp

    if (present(acache)) then
      aptr => acache
    else
      allocate (local_acache)
      if (ref%nat .ne. mol%nat) then
        error stop 'Unequal molecule size in compute_linear_sum_assignment()'
      end if
      nat = max(ref%nat,mol%nat)
      call local_acache%allocate(nat,nat,.true.)
      aptr => local_acache
    end if

    !> Compute the cost matrix, which is simply the distance matrix
    !> between the two molecules.
    !> To avoid computational overhead we can skip the square root.
    !> It won't affect the result
    !> Also, since aptr%Cost is a flattened matrix, we only fill
    !> the first rnknat**2 entries
    rnknat = ngroups(targetrank)
    ii = 0
    do i = 1,ref%nat
      if (ranks(i,1) .ne. targetrank) cycle
      ii = ii+1
      iwork2(ii,1) = i !> mapping using the first column of iwork2
      jj = 0
      do j = 1,mol%nat
        if (ranks(j,2) .ne. targetrank) cycle
        jj = jj+1
        dists(:) = (ref%xyz(:,i)-mol%xyz(:,j))**2 !> use i and j
        aptr%Cost(jj+(ii-1)*rnknat) = sum(dists)
      end do
    end do

    if (debug) then
      write (*,*) 'target rank',targetrank,'# atoms',rnknat
    end if

    call lsap(aptr,rnknat,rnknat,.false.,iostatus)

    !> paasing back the determined order as second column of iwork2
    if (iostatus == 0) then
      if (debug) then
        do i = 1,rnknat
          write (*,*) iwork2(aptr%a(i),1),'-->',iwork2(aptr%b(i),1)
        end do
      end if
      do i = 1,rnknat
        jj = aptr%a(i)
        ii = aptr%b(i)
        if(ii == -1 .or. jj == -1) cycle  !> cycle bad assignments
        val0 = val0+aptr%Cost(jj+(ii-1)*rnknat)
        iwork2(i,2) = iwork2(aptr%b(i),1)
      end do
    else
      !> in the unlikely case we have a failure of the LSAP
      !> we do just a 1:1 mapping, just so that the algo doesn't crash
      iwork2(1:rnknat,2) = iwork2(1:rnknat,1)
    end if

  end subroutine compute_linear_sum_assignment

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

  function checkranks(nat,ranks1,ranks2) result(yesno)
!***********************************************************************
!* Check two rank arrays to see if we have the same amount of
!* atoms in the same ranks (a condition to bein able to work with them)
!***********************************************************************
    implicit none
    logical :: yesno
    integer,intent(in) :: nat
    integer,intent(in) :: ranks1(nat)
    integer,intent(in) :: ranks2(nat)
    integer :: ii,jj,maxrank1,maxrank2
    integer :: count1,count2
    yesno = .false.

    maxrank1 = maxval(ranks1)
    maxrank2 = maxval(ranks2)
    !> different maxranks, so we can't have the same and return
    if (maxrank1 .ne. maxrank2) return

    do ii = 1,maxrank1
      count1 = 0
      count2 = 0
      do jj = 1,nat
        if (ranks1(jj) .eq. ii) count1 = count1+1
        if (ranks2(jj) .eq. ii) count2 = count2+1
      end do
      !> not the same amount of atoms in rank ii, return from function
      if (count1 .ne. count2) return
    end do

    !> if we reach this point we can assume the given ranks are o.k.
    yesno = .true.
  end function checkranks

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
