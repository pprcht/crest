module hessian_reconstruct
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public cashed_hessian

  type :: cashed_hessian

    integer :: steps = 10
    real(wp),allocatable :: gradient(:,:,:)
    real(wp),allocatable :: coords(:,:,:)
    real(wp),allocatable :: energy(:)
    real(wp),allocatable :: s(:,:),y(:,:),B(:,:),H(:,:),p(:),rho(:),V(:,:,:),I(:,:)
    integer,allocatable :: order(:),natm
    integer :: stepcount = 0

  contains

    procedure :: alloc => cashed_hessian_allocate
    procedure :: dealloc => cashed_hessian_deallocate
    procedure :: update => update_cashed_hessian
    procedure :: construct_hessian_lbfgs
    procedure :: compute_intermediates

  end type cashed_hessian

contains

  subroutine cashed_hessian_allocate(self,N,steps)
    integer,intent(in) :: N,steps
    class(cashed_hessian),intent(inout) :: self

    self%steps = steps
    allocate (self%gradient(steps,3,N))
    allocate (self%coords(steps,3,N))
    allocate (self%energy(steps))
    allocate (self%order(steps))
    allocate (self%s(self%steps-1,3*self%natm))
    allocate (self%y(self%steps-1,3*self%natm))
    allocate (self%p(self%steps-1))
    allocate (self%rho(self%steps-1))
    allocate (self%V(self%steps-1,3*self%natm,3*self%natm))
    allocate (self%I(3*self%natm,3*self%natm))
    self%natm = N

    self%order(:) = 0

  end subroutine cashed_hessian_allocate

  subroutine cashed_hessian_deallocate(self)
    class(cashed_hessian),intent(inout) :: self

    if (allocated(self%gradient)) deallocate (self%gradient)
    if (allocated(self%coords)) deallocate (self%coords)
    if (allocated(self%energy)) deallocate (self%energy)
    if (allocated(self%order)) deallocate (self%order)

  end subroutine cashed_hessian_deallocate

  subroutine update_cashed_hessian(self,gradient,energy,coords)
    class(cashed_hessian),intent(inout) :: self
    real(wp),intent(in) :: gradient(:,:),energy,coords(:,:)
    integer :: idx

    self%stepcount = self%stepcount+1
    idx = minloc(self%order,1)
    self%order(idx) = self%stepcount
    self%gradient(idx,:,:) = gradient
    self%energy(idx) = energy
    self%coords(idx,:,:) = coords

  end subroutine update_cashed_hessian

  recursive subroutine construct_hessian_lbfgs(self,n)
    class(cashed_hessian),intent(inout) :: self
    integer,intent(in) :: n
    real(wp),allocatable :: temp(:,:)

    allocate (temp(3*self%natm,3*self%natm))
    if (n == 0) then
      call self%compute_intermediates()
      allocate (self%B(3*self%natm,3*self%natm))
      self%B = self%I
    else
      call self%construct_hessian_lbfgs(n-1)
      temp = matmul(matmul(TRANSPOSE(self%V(n,:,:)),self%B),self%V(n,:,:))-self%p(n)*(matmul(reshape(self%s(n,:), [3*self%natm,1]),reshape(self%s(n,:), [1,3*self%natm])))
      self%B = temp
    end if

  end subroutine construct_hessian_lbfgs

  subroutine compute_intermediates(self)
    class(cashed_hessian),intent(inout) :: self
    integer :: i,j,k
    real(wp), allocatable :: tmp(:),tmp_coords(:,:),tmp_grads(:,:)

    allocate (tmp_coords(self%steps,3*self%natm))
    allocate (tmp_grads(self%steps,3*self%natm))
    allocate (tmp(self%steps))
    allocate (self%s(self%steps-1,3*self%natm))
    allocate (self%y(self%steps-1,3*self%natm))
    allocate (self%p(self%steps-1))
    allocate (self%rho(self%steps-1))
    allocate (self%V(self%steps-1,3*self%natm,3*self%natm))
    allocate (self%I(3*self%natm,3*self%natm))

    tmp = self%order
    self%I = 0.0_wp

    do k = 1,3*self%natm
      self%I(k,k) = 1.0_wp
    end do

    tmp_coords = reshape(self%coords,[self%steps,3*self%natm])
    tmp_grads = reshape(self%gradient,[self%steps,3*self%natm])

    if (minval(tmp) == 0) then
      print*,"ERROR: Number of recursive steps for hessian reconstruction larger than number of geoemtry optimization steps!"
    else
      do i = 1,self%steps
        if (i == 1) then
          j = minloc(tmp,1)
          tmp(j) = HUGE(tmp(j))
        else
          j = minloc(tmp,1)
          if (j == 1) then
            self%s(i-1,:) = tmp_coords(j,:)-tmp_coords(self%steps,:)
            self%y(i-1,:) = tmp_grads(j,:)-tmp_grads(self%steps,:)
          else
            self%s(i-1,:) = tmp_coords(j,:)-tmp_coords(j-1,:)
            self%y(i-1,:) = tmp_grads(j,:)-tmp_grads(j-1,:)
          end if
          self%p(i-1) = 1/(dot_product(self%y(i-1,:),self%s(i-1,:)))
          self%V(i-1,:,:) = self%I-self%p(i-1)*(matmul(reshape(self%y(i-1,:), [3*self%natm,1]),reshape(self%s(i-1,:), [1,3*self%natm])))
          tmp(j) = HUGE(tmp(j))
        end if
      end do
    end if

  end subroutine compute_intermediates

end module hessian_reconstruct
