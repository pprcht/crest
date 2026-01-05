module construct_list
  !**************************************************
  !* Bookkeeping module for reconstructing moleucles
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

  type :: construct_layer
    integer :: nnodes = 0
    type(coord),allocatable :: node(:)
    !> initial/reconstructed molecule
    type(coord) :: mol
  end type construct_layer

  type :: construct_heap
    integer :: nlayer = 0
    type(construct_layer),allocatable :: layer(:)

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

end module construct_list
