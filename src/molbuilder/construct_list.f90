module construct_list
  !**************************************************
  !* Bookkeeping module for reconstructing molecules
  !**************************************************
  use crest_parameters
  use strucrd
  use quicksort_interface,only:qqsorti
  implicit none
  private

  type :: split_atms
    integer :: natms = 0
    integer,allocatable :: atms(:)
  end type split_atms

  type :: origin_map
    integer :: natms = 0
    integer,allocatable :: map(:)
  end type origin_map

  type :: construct_layer
    integer :: id = 0
    integer :: parent = 0
    integer :: parentnode = 0
    integer :: nnodes = 0
    type(coord),allocatable :: node(:)
    integer,allocatable :: childlayer(:)
    integer,allocatable :: alignmap(:,:)
    integer,allocatable :: ncapped(:)
    integer,allocatable :: position_mapping(:,:)
    type(origin_map),allocatable :: origin(:)
    !> initial/reconstructed molecule
    type(coord) :: mol
  end type construct_layer

  type :: construct_queue
    integer :: layer = 0
    integer :: node = 0
    integer :: duplicate_of_queue = 0
    character(len=:),allocatable :: workdir
    character(len=:),allocatable :: file
  end type construct_queue

  type :: construct_heap
    !> The gerneral mapping
    integer :: nlayer = 0
    type(construct_layer),allocatable :: layer(:)
    !> the queue to treat endpoints
    integer :: nqueue = 0
    type(construct_queue),allocatable :: queue(:)
    !> originial directory
    character(len=:),allocatable :: origindir
  contains
    procedure map_origins_for_layer
    procedure find_current_position
    procedure count_endpoints
    procedure setup_queue
  end type construct_heap

  !> exported types
  public :: split_atms
  public :: construct_heap
  !> exported helper functions
  public :: add_to_splitqueue

!=============================================================================!
contains  !> MODULE PROCEDURES START HERE
!=============================================================================!

  subroutine add_to_splitqueue(splitqueue,raw_split)
    implicit none
    type(split_atms),intent(inout),allocatable :: splitqueue(:)
    integer,intent(in) :: raw_split(:)
    type(split_atms) :: tmpsplit
    type(split_atms),allocatable :: tmpsplitqueue(:)
    integer :: ii,n,nall,nallnew
    logical :: duplicate

    n = size(raw_split,1)
    allocate (tmpsplit%atms(n))
    tmpsplit%natms = n
    tmpsplit%atms(:) = raw_split(:)
    call qqsorti(tmpsplit%atms,1,n)

    if (.not.allocated(splitqueue)) then
      allocate (splitqueue(1))
      splitqueue(1) = tmpsplit
    else
      nall = size(splitqueue,1)
      nallnew = nall+1
      duplicate = .false.
      do ii = 1,nall
        if (all(splitqueue(ii)%atms(:) .eq. tmpsplit%atms(:))) then
          duplicate = .true.
          exit
        end if
      end do
      if (.not.duplicate) then
        allocate (tmpsplitqueue(nallnew))
        do ii = 1,nall
          tmpsplitqueue(ii) = splitqueue(ii)
        end do
        tmpsplitqueue(nallnew) = tmpsplit
        call move_alloc(tmpsplitqueue,splitqueue)
      end if
    end if

  end subroutine add_to_splitqueue

  recursive function find_original_atom(atom,heap,targetlayer,targetnode) result(this)
    implicit none
    type(construct_heap),intent(in) :: heap
    integer,intent(in) :: targetlayer,targetnode,atom
    integer :: this
    integer :: ii,jj,kk

    this = 0
    kk = 0
    do ii = 1,size(heap%layer(targetlayer)%position_mapping,1)
      if (atom == heap%layer(targetlayer)%position_mapping(ii,targetnode)) then
        kk = ii
        exit
      end if
    end do
    if (heap%layer(targetlayer)%parent == 0.or.kk == 0) then
      this = kk
    else
      ii = heap%layer(targetlayer)%parent
      jj = heap%layer(targetlayer)%parentnode
      this = find_original_atom(kk,heap,ii,jj)
    end if
  end function find_original_atom

  subroutine map_origins_for_layer(heap,targetlayer)
    implicit none
    class(construct_heap),intent(inout) :: heap
    integer,intent(in) :: targetlayer
    integer :: ii,jj,kk,nat
    logical,parameter :: debug = .false.
    if (targetlayer < 1.or.targetlayer > heap%nlayer) return
    if (.not.allocated(heap%layer(targetlayer)%node)) return
    associate (layer => heap%layer(targetlayer))
      if (allocated(layer%origin)) deallocate (layer%origin)
      allocate (layer%origin(layer%nnodes))
      do ii = 1,layer%nnodes
        nat = layer%node(ii)%nat
        allocate (layer%origin(ii)%map(nat),source=0)
        layer%origin(ii)%natms = nat
        do jj = 1,nat
          layer%origin(ii)%map(jj) = find_original_atom(jj,heap,targetlayer,ii)
        end do
        if (debug) then
          write (stdout,*) 'Original atom positions for fragment',ii,'of layer',targetlayer
          do jj = 1,nat
            write (stdout,*) 'frag.atm.',jj,'<-- original:',layer%origin(ii)%map(jj)
          end do
        end if
      end do
    end associate
  end subroutine map_origins_for_layer

  subroutine find_current_position(heap,atom,targetlayer,targetnode,current)
    implicit none
    class(construct_heap),intent(in) :: heap
    integer,intent(in) :: atom,targetlayer,targetnode
    integer,intent(out) :: current
    integer :: ii,jj
    current = 0
    associate (layer => heap%layer(targetlayer))
      do ii = 1,layer%node(targetnode)%nat
        jj = layer%origin(targetnode)%map(ii)
        if (jj == atom) then
          current = ii
          exit
        end if
      end do
    end associate
  end subroutine find_current_position

  function count_endpoints(heap) result(nendpoints)
    implicit none
    class(construct_heap) :: heap
    integer :: nendpoints
    integer :: ii,jj,kk
    nendpoints = 0
    do ii = 1,heap%nlayer
      if (.not.allocated(heap%layer(ii)%childlayer)) then
        nendpoints = nendpoints+heap%layer(ii)%nnodes
      else
        do jj = 1,heap%layer(ii)%nnodes
          if (heap%layer(ii)%childlayer(jj) == 0) then
            nendpoints = nendpoints+1
          end if
        end do
      end if
    end do
  end function count_endpoints

  subroutine setup_queue(heap)
    implicit none
    class(construct_heap) :: heap
    integer :: nqueue,kk,ii,jj

    nqueue = heap%count_endpoints()
    heap%nqueue = nqueue
    allocate (heap%queue(nqueue))
    
    kk = 0
    do ii = 1,heap%nlayer
      if (.not.allocated(heap%layer(ii)%childlayer)) then
        do jj = 1,heap%layer(ii)%nnodes
          kk = kk+1
          heap%queue(kk)%layer = ii
          heap%queue(kk)%node = jj
        end do
      else
        do jj = 1,heap%layer(ii)%nnodes
          if (heap%layer(ii)%childlayer(jj) == 0) then
            kk = kk+1
            heap%queue(kk)%layer = ii
            heap%queue(kk)%node = jj
          end if
        end do
      end if
    end do

  end subroutine setup_queue

end module construct_list
