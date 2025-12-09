module hessian_reconstruct
  use iso_fortran_env,only:wp => real64
  use hessupdate_module
  use optimize_maths
  implicit none
  private

  public cashed_hessian,invert_matrix

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

  contains

    procedure :: alloc => cashed_hessian_allocate
    procedure :: dealloc => cashed_hessian_deallocate
    procedure :: update => update_cashed_hessian
    procedure :: construct_hessian_lbfgs
    procedure :: compute_intermediates
    procedure :: construct_hessian_bfgs

  end type cashed_hessian

contains

  subroutine cashed_hessian_allocate(self,N,steps,hguess) !> maybe make keywords optional later
    integer,intent(in) :: N,steps
    class(cashed_hessian),intent(inout) :: self
    real(wp),intent(in) :: hguess

    self%steps = steps
    self%hguess = hguess
    self%natm = N
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
    !print*, coords
    self%stepcount = self%stepcount+1
    !print*, self%order
    !print*, coords(:,:)
    idx = minloc(self%order,1)
    self%order(idx) = self%stepcount
    self%gradient(idx,:,:) = gradient
    self%energy(idx) = energy
    self%coords(idx,:,:) = coords
    !if (idx==1) then
    !  print*, self%coords(1,:,:)
    !  print*, self%coords(2,:,:)
    !  print*, self%coords(3,:,:)
    !endif
    !PRINT*, self%order
    !print*, self%coords(1,:,:)
    !print*, self%coords(2,:,:)
    print*,size(self%coords,dim=1)

  end subroutine update_cashed_hessian

  subroutine construct_hessian_bfgs(self)
    class(cashed_hessian),intent(inout) :: self
    integer :: i,j,k,nat3
    real(wp),allocatable :: tmp(:),tmp_coords(:,:),tmp_grads(:,:),hess(:),dx(:)
    real(wp) :: gnorm

    nat3 = 3*self%natm

    allocate (tmp_coords(self%steps,nat3))
    allocate (tmp_grads(self%steps,nat3))
    allocate (tmp(self%steps))
    allocate (hess(nat3*(nat3+1)/2))
    allocate (dx(nat3))

    tmp = self%order

    tmp_coords = reshape(self%coords, [self%steps,nat3])
    tmp_grads = reshape(self%gradient, [self%steps,nat3])

    if (minval(tmp) == 0) then
      print*,"ERROR: Number of recursive steps for hessian reconstruction larger than number of geoemtry optimization steps!"
    else
      do i = 1,self%steps
        if (i == 1) then
          j = minloc(tmp,1)
          tmp(j) = HUGE(tmp(j))
          do k = 1,nat3
            self%hguess_mat(k,k) = self%hguess
          end do
          call dsqtoh(nat3,self%hguess_mat,hess)
        else
          j = minloc(tmp,1)
          if (j == 1) then
            dx = tmp_coords(j,:)-tmp_coords(self%steps,:)
            call bfgs(nat3,gnorm,tmp_grads(j,:),tmp_grads(self%steps,:),dx,hess)
          else
            dx = tmp_coords(j,:)-tmp_coords(j-1,:)
            call bfgs(nat3,gnorm,tmp_grads(j,:),tmp_grads(j-1,:),dx,hess)
          end if
          tmp(j) = HUGE(tmp(j))
        end if
      end do
    end if

    call dhtosq(nat3,self%B,hess)

  end subroutine construct_hessian_bfgs

  recursive subroutine construct_hessian_lbfgs(self,n)  !> refactor this to reduce memory by
    class(cashed_hessian),intent(inout) :: self         !> computing intermediates within this routine
    integer,intent(in) :: n
    real(wp),allocatable :: temp(:,:)
    !real(wp), allocatable :: test_mat(:,:)

    !allocate (test_mat(3*self%natm,3*self%natm))

    allocate (temp(3*self%natm,3*self%natm))
    if (n == 0) then
      call self%compute_intermediates()
      allocate (self%B(3*self%natm,3*self%natm))
      self%B = self%hguess
    else
      call self%construct_hessian_lbfgs(n-1)
      temp = matmul(matmul(TRANSPOSE(self%V(n,:,:)),self%B),self%V(n,:,:))+self%p(n)*(matmul(reshape(self%s(n,:), [3*self%natm,1]),reshape(self%s(n,:), [1,3*self%natm])))
      self%B = temp
      print*
      print*,"updated Hessian number",N
      print*
      print*,temp(1,:)
    end if

  end subroutine construct_hessian_lbfgs

  subroutine compute_intermediates(self)
    class(cashed_hessian),intent(inout) :: self
    integer :: i,j,k,l
    real(wp),allocatable :: tmp(:),tmp_coords(:,:),tmp_grads(:,:)
    real(wp),allocatable :: temp_mat(:,:)

    allocate (temp_mat(3*self%natm,3*self%natm))

    allocate (tmp_coords(self%steps,3*self%natm))
    allocate (tmp_grads(self%steps,3*self%natm))
    allocate (tmp(self%steps))

    tmp = self%order
    self%I = 0.0_wp

    do k = 1,3*self%natm
      self%I(k,k) = 1.0_wp
    end do

    self%hguess_mat = 0.0_wp

    do l = 1,3*self%natm
      self%hguess_mat = self%hguess
    end do

    tmp_coords = reshape(self%coords, [self%steps,3*self%natm])
    tmp_grads = reshape(self%gradient, [self%steps,3*self%natm])

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
          self%V(i-1,:,:) = (self%I(:,:))-(self%p(i-1)*(matmul(reshape(self%y(i-1,:), [3*self%natm,1]),reshape(self%s(i-1,:), [1,3*self%natm]))))
          tmp(j) = HUGE(tmp(j))
          !temp_mat = self%p(i-1)*(matmul(reshape(self%y(i-1,:), [3*self%natm,1]),reshape(self%s(i-1,:), [1,3*self%natm])))
        end if
      end do
    end if

  end subroutine compute_intermediates

  ! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
  function invert_matrix(A) result(Ainv)
    real(wp),dimension(:,:),intent(in) :: A
    real(wp),dimension(size(A,1),size(A,2)) :: Ainv

    real(wp),dimension(size(A,1)) :: work  ! work array for LAPACK
    integer,dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n,info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n,n,Ainv,n,ipiv,info)

    if (info /= 0) then
      stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n,Ainv,n,ipiv,work,n,info)

    if (info /= 0) then
      stop 'Matrix inversion failed!'
    end if
  end function invert_matrix

end module hessian_reconstruct
