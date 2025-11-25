module hessian_reconstruct
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  public cashed_hessian

  type :: cashed_hessian

  integer :: steps = 10
  real(wp), allocatable :: gradient(:,:,:)
  real(wp), allocatable :: coords(:,:,:)
  real(wp), allocatable :: energy(:)
  real(wp), allocatable :: order(:)
  integer :: stepcount = 0

  contains

  procedure :: alloc => cashed_hessian_allocate
  procedure :: dealloc => cashed_hessian_deallocate
  procedure :: update => update_cashed_hessian

  end type cashed_hessian 

contains

subroutine cashed_hessian_allocate(self, N, steps)
integer, intent(in) :: N, steps
class(cashed_hessian) :: self

self%steps = steps
allocate(self%gradient(steps,3,N))
allocate(self%coords(steps,3,N))
allocate(self%energy(steps))
allocate(self%order(steps))

self%order(:) = 0.0_wp

end subroutine cashed_hessian_allocate

subroutine cashed_hessian_deallocate(self)
class(cashed_hessian) :: self

if(allocated(self%gradient)) deallocate(self%gradient)
if(allocated(self%coords)) deallocate(self%coords)
if(allocated(self%energy)) deallocate(self%energy)
if(allocated(self%order)) deallocate(self%order)

end subroutine cashed_hessian_deallocate

subroutine update_cashed_hessian(self, gradient, energy, coords)
class(cashed_hessian) :: self
real(wp), intent(in) :: gradient(:,:), energy, coords(:,:)
integer :: idx

self%stepcount = self%stepcount + 1
idx = minloc(self%order,1)
self%order(idx) = self%stepcount
self%gradient(idx,:,:) = gradient
self%energy(idx) = energy
self%coords(idx,:,:) = coords

end subroutine update_cashed_hessian

end module hessian_reconstruct