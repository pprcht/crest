module hessian_reconstruct
  use iso_fortran_env,only:wp => real64
  use hessupdate_module
  use optimize_maths
  use crest_parameters
  implicit none
  private

  public cashed_hessian

  type :: cashed_hessian

    integer :: steps = 10
    real(wp),allocatable :: gradient(:,:,:)
    real(wp),allocatable :: coords(:,:,:)
    real(wp),allocatable :: energy(:)
    real(wp),allocatable :: s(:,:),y(:,:),B(:,:),H(:,:),Hinv(:,:),p(:),rho(:),V(:,:,:),I(:,:)
    integer,allocatable :: order(:),natm
    integer :: stepcount = 0
    real(wp) :: hguess = 0.02_wp
    real(wp),allocatable ::hguess_mat(:,:)
    logical :: track_step = .true.
    integer :: initialize_type = 0

  contains

    procedure :: alloc => cashed_hessian_allocate
    procedure :: dealloc => cashed_hessian_deallocate
    procedure :: update => update_cashed_hessian
    procedure :: construct_hessian_bfgs

  end type cashed_hessian

contains

  subroutine cashed_hessian_allocate(self,N,steps,hguess,initialize_type) !> maybe make keywords optional later
    integer,intent(in) :: N,steps, initialize_type
    class(cashed_hessian),intent(inout) :: self
    real(wp),intent(in) :: hguess

    self%steps = steps
    self%hguess = hguess
    self%natm = N
    self%initialize_type = initialize_type
    allocate (self%gradient(steps,3,N))
    allocate (self%coords(steps,3,N))
    allocate (self%energy(steps))
    allocate (self%order(steps))
    allocate (self%s(self%steps-1,3*N))
    allocate (self%y(self%steps-1,3*N))
    allocate (self%p(self%steps-1))
    allocate (self%rho(self%steps-1))
    allocate (self%V(self%steps-1,3*N,3*N))
    allocate (self%I(3*N,3*N))
    allocate (self%hguess_mat(3*N,3*N))
    allocate (self%H(3*N,3*N))
    allocate (self%Hinv(3*N,3*N))
    allocate (self%B(3*N,3*N))

    self%order(:) = 0

  end subroutine cashed_hessian_allocate

  subroutine cashed_hessian_deallocate(self)
    class(cashed_hessian),intent(inout) :: self

    if (allocated(self%gradient)) deallocate (self%gradient)
    if (allocated(self%coords)) deallocate (self%coords)
    if (allocated(self%energy)) deallocate (self%energy)
    if (allocated(self%order)) deallocate (self%order)
    if (allocated(self%s)) deallocate (self%s)
    if (allocated(self%y)) deallocate (self%y)
    if (allocated(self%p)) deallocate (self%p)
    if (allocated(self%rho)) deallocate (self%rho)
    if (allocated(self%V)) deallocate (self%V)
    if (allocated(self%I)) deallocate (self%I)

  end subroutine cashed_hessian_deallocate

  subroutine update_cashed_hessian(self,gradient,energy,coords)
    class(cashed_hessian),intent(inout) :: self
    real(wp),intent(in) :: gradient(:,:),energy,coords(:,:)
    integer :: idx,i

    self%stepcount = self%stepcount+1
    idx = minloc(self%order,1)
    self%order(idx) = self%stepcount
    self%gradient(idx,:,:) = gradient
    self%energy(idx) = energy
    self%coords(idx,:,:) = coords

  end subroutine update_cashed_hessian

  subroutine construct_hessian_bfgs(self)
    class(cashed_hessian),intent(inout) :: self
    integer :: i,j,k,nat3
    real(wp),allocatable :: tmp(:),tmp_coords(:,:),tmp_grads(:,:),hess(:),dx(:)
    real(wp) :: gnorm
    integer :: unit,iter,made_iters

    nat3 = 3*self%natm

    allocate (tmp_coords(self%steps,nat3))
    allocate (tmp_grads(self%steps,nat3))
    allocate (tmp(self%steps))
    allocate (hess(nat3*(nat3+1)/2))
    allocate (dx(nat3))

    tmp = self%order

    tmp_coords = reshape(self%coords, [self%steps,nat3])
    tmp_grads = reshape(self%gradient, [self%steps,nat3])

    made_iters = self%steps

    call dsqtoh(nat3,self%hguess_mat,hess) !> Here, Hessian is packed and transferred to hess

    if (minval(tmp) == 0) then !> Implement keyword like exact HU that kills the process
      made_iters = maxval(tmp) !> if made_iters<steps
      write (stdout,*) "Requsted Number of reconstruction steps is",self%steps, &
      & "but only",made_iters,"geometry optimization steps were made!"
      write (stdout,*) "Hessian is reconstructed with",made_iters,"update steps only!"

      do while (minval(tmp) == 0)
        j = minloc(tmp,1)
        tmp(j) = HUGE(tmp(j))
      end do
    end if

    do i = 1,made_iters
      if (i == 1) then
        j = minloc(tmp,1)
        tmp(j) = HUGE(tmp(j))
      else
        j = minloc(tmp,1) !> This only happens if made_iters>steps
        if (j == 1) then  !> => Not affected if too many steps requested
          dx = tmp_coords(j,:)-tmp_coords(self%steps,:)
          call bfgs(nat3,gnorm,tmp_grads(j,:),tmp_grads(self%steps,:),dx,hess)
        else
          dx = tmp_coords(j,:)-tmp_coords(j-1,:)
          call bfgs(nat3,gnorm,tmp_grads(j,:),tmp_grads(j-1,:),dx,hess)
        end if
        tmp(j) = HUGE(tmp(j))
      end if
    end do

    call dhtosq(nat3,self%B,hess)

  end subroutine construct_hessian_bfgs

end module hessian_reconstruct
