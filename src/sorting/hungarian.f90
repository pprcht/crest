module hungarian_module
!************************************************************
!* Implementations of the Hungarian (Kuhn-Munkres) Algorithm
!* in O(n³) time (Edmons & Karp / Tomizawa).
!*
!* Implemented in single precision with a cache to
!* circumvent repeated memory allocation.
!*
!* Also includes some wrappers for standalone use
!************************************************************
  use iso_fortran_env,sp => real32,wp => real64
  implicit none
  private

  public :: hungarian
  interface hungarian
    module procedure hungarian_cached
    module procedure hungarian_wrap_int
    module procedure hungarian_wrap_sp
    module procedure hungarian_wrap_wp
  end interface hungarian

  real(sp),parameter,private :: inf = huge(1.0_sp) !> Use huge intrinsic for large numbers
  integer,parameter,private  :: infi = huge(1) !> Use huge intrinsic for large numbers

  public :: assignment_cache
  type :: assignment_cache
    integer :: J,W
    real(sp),allocatable :: Cost(:,:)  !> Cost(J,W)
    real(sp),allocatable :: answers(:) !> answers(J)
    integer,allocatable  :: job(:)     !> job(W+1)
    !> Workspace
    real(sp),allocatable  :: ys(:)     !> ys(J)
    real(sp),allocatable  :: yt(:)     !> yt(W+1)
    real(sp),allocatable  :: Ct(:,:)   !> Ct(W,J)
    real(sp),allocatable  :: min_to(:) !> min_to(W+1)
    integer,allocatable  :: prv(:)     !> prv(W+1)
    logical,allocatable  :: in_Z(:)    !> in_Z(W+1)
  contains
    procedure :: allocate => allocate_assignment_cache
    procedure :: deallocate => deallocate_assignment_cache
  end type assignment_cache

  interface ckmin
    module procedure ckmin_int
    module procedure ckmin_sp
  end interface ckmin

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine allocate_assignment_cache(self,J,W)
    implicit none
    class(assignment_cache),intent(inout) :: self
    integer,intent(in) :: J,W
    !> Store dimensions
    self%J = J
    self%W = W
    !> Allocate arrays based on input dimensions
    allocate (self%Cost(J,W))
    allocate (self%answers(J))
    allocate (self%job(W+1))
    !> Allocate workspace arrays
    allocate (self%ys(J))
    allocate (self%yt(W+1))
    allocate (self%Ct(W,J))
    allocate (self%min_to(W+1))
    allocate (self%prv(W+1))
    allocate (self%in_Z(W+1))
  end subroutine allocate_assignment_cache

  subroutine deallocate_assignment_cache(self)
    implicit none
    class(assignment_cache),intent(inout) :: self
    ! Deallocate arrays if they are allocated
    if (allocated(self%Cost)) deallocate (self%Cost)
    if (allocated(self%answers)) deallocate (self%answers)
    if (allocated(self%job)) deallocate (self%job)
    if (allocated(self%ys)) deallocate (self%ys)
    if (allocated(self%yt)) deallocate (self%yt)
    if (allocated(self%Ct)) deallocate (self%Ct)
    if (allocated(self%min_to)) deallocate (self%min_to)
    if (allocated(self%prv)) deallocate (self%prv)
    if (allocated(self%in_Z)) deallocate (self%in_Z)
  end subroutine deallocate_assignment_cache

!========================================================================================!

  logical function ckmin_int(a,b) result(yesno)
    !> Helper function to compute the minimum and update
    integer,intent(inout) :: a
    integer,intent(in) :: b
    yesno = .false.
    if (b < a) then
      a = b
      yesno = .true.
    end if
  end function ckmin_int

  logical function ckmin_sp(a,b) result(yesno)
    !> Helper function to compute the minimum and update
    real(sp),intent(inout) :: a
    real(sp),intent(in) :: b
    yesno = .false.
    if (b < a) then
      a = b
      yesno = .true.
    end if
  end function ckmin_sp

  subroutine hungarian_cached(cache,J,W)
    !****************************************************************
    !* Hungarian algorithm implementation to solve an assignment
    !* problem in O(n³) time.
    !* This implementation refers to a cache, which is created
    !* to avoid repeated memory allocation.
    !* Passing J and W explicitly enables reuse of memory
    !* for smaller sub-problems (i.e. cache%J >= J, W accoridingly)
    !*
    !* Inputs (all within cache, except J and W):
    !*   C(J, W) - Cost matrix of dimensions J-by-W,
    !*             where C(jj, ww) is the cost to assign
    !*             jj-th job to ww-th worker
    !*   J       - Number of jobs
    !*   W       - Number of workers
    !* Outputs (all within cache):
    !*   answers(J) - Vector of length J, where answers(jj) is
    !*                the minimum cost to assign the first jj
    !*                jobs to distinct workers
    !*   job(W+1)   - Vector where job(ww) is the job assigned to
    !*                the ww-th worker (or -1 if no job is assigned)
    !****************************************************************
    integer,intent(in) :: J,W
    type(assignment_cache),intent(inout) :: cache
    integer  :: jj_cur,ww_cur,jj,ww_next,ww
    real(sp) :: delta

    !> IMPORTANT: associate to have shorter variable names
    associate (C => cache%Cost, &
              & answers => cache%answers, &
              & job => cache%job, &
              & ys => cache%ys, &
              & yt => cache%yt, &
              & Ct => cache%Ct, &
              & min_to => cache%min_to, &
              & prv => cache%prv, &
              & in_Z => cache%in_Z)

      job = -1
      ys = 0
      yt = 0
      Ct = transpose(C)

      do jj_cur = 1,J   !> O(n¹)
        ww_cur = W+1
        job(ww_cur) = jj_cur
        min_to = inf
        prv = -1
        in_Z = .false.

        do while (job(ww_cur) /= -1)  !> O(n¹) -> O(n²)
          in_Z(ww_cur) = .true.
          jj = job(ww_cur)
          delta = inf
          do ww = 1,W !> O(n²) -> O(n³)
            if (.not.in_Z(ww)) then
              if (ckmin(min_to(ww),Ct(ww,jj)-ys(jj)-yt(ww))) then
                prv(ww) = ww_cur
              end if
              if (ckmin(delta,min_to(ww))) then
                ww_next = ww
              end if
            end if
          end do

          do ww = 1,W+1
            if (in_Z(ww)) then
              ys(job(ww)) = ys(job(ww))+delta
              yt(ww) = yt(ww)-delta
            else
              min_to(ww) = min_to(ww)-delta
            end if
          end do
          ww_cur = ww_next
        end do

        !> Update assignments along alternating path
        do while (ww_cur /= W+1)
          job(ww_cur) = job(prv(ww_cur))
          ww_cur = prv(ww_cur)
        end do

        answers(jj_cur) = -yt(W+1)
      end do

    end associate
  end subroutine hungarian_cached

!========================================================================================!

  subroutine hungarian_wrap_int(C,J,W,answers,job)
    !*********************************************
    !* Wrapper for integer precision
    !*********************************************
    implicit none
    integer,intent(in) :: C(J,W)
    integer,intent(in) :: J,W
    integer,intent(out) :: answers(J)
    integer,intent(out) :: job(W+1)
    type(assignment_cache) :: cache

    call cache%allocate(J,W)
    cache%Cost(1:J,1:W) = real(C(1:J,1:W),sp)
    call hungarian_cached(cache,J,W)

    answers(1:J) = nint(cache%answers(1:J))
    job(1:W+1) = cache%job(1:W+1)
    call cache%deallocate()
  end subroutine hungarian_wrap_int

  subroutine hungarian_wrap_sp(C,J,W,answers,job)
    !*********************************************
    !* Wrapper for single precision
    !*********************************************
    implicit none
    real(sp),intent(in) :: C(J,W)
    integer,intent(in) :: J,W
    real(sp),intent(out) :: answers(J)
    integer,intent(out) :: job(W+1)
    type(assignment_cache) :: cache

    call cache%allocate(J,W)
    cache%Cost(1:J,1:W) = C(1:J,1:W)
    call hungarian_cached(cache,J,W)

    answers(1:J) = cache%answers(1:J)
    job(1:W+1) = cache%job(1:W+1)
    call cache%deallocate()
  end subroutine hungarian_wrap_sp

  subroutine hungarian_wrap_wp(C,J,W,answers,job)
    !*********************************************
    !* Wrapper for double precision
    !*********************************************
    implicit none
    real(wp),intent(in) :: C(J,W)
    integer,intent(in) :: J,W
    real(wp),intent(out) :: answers(J)
    integer,intent(out) :: job(W+1)
    type(assignment_cache) :: cache

    call cache%allocate(J,W)
    cache%Cost(1:J,1:W) = real(C(1:J,1:W),sp)
    call hungarian_cached(cache,J,W)

    answers(1:J) = real(cache%answers(1:J),wp)
    job(1:W+1) = cache%job(1:W+1)
    call cache%deallocate()
  end subroutine hungarian_wrap_wp

!========================================================================================!
!========================================================================================!
end module hungarian_module
