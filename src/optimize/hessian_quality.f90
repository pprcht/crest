!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht, Lukas Rindt
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!> Practical diagnostics tools for assessing and improving approximate Hessians
!> from geometry optimisations, targeted at thermochemistry applications.

module hessian_quality
  use crest_parameters
  use strucrd
  use crest_calculator
  implicit none
  private
  public :: step_monitor,monitor_init,monitor_add_step,monitor_should_freeze, &
            mode_quality_analysis,selective_hessian_repair

  !>────────────────────────────────────────────────────────────────────────────
  !> Incremental step monitor.
  !>
  !> Stores the orthonormal basis Q of all accepted steps (those with
  !> sufficient information gain).  Operates as a "growing QR"; at each
  !> step we orthogonalise against Q and either accept (expand Q by one
  !> column) or reject (information gain below threshold).
  type :: step_monitor
    integer  :: n3 = 0              !> flattened dimension (3 * nat)
    integer  :: max_cols = 0        !> allocated capacity
    integer  :: rank = 0            !> current rank  (# accepted directions)
    integer  :: n_total = 0         !> total steps seen (including rejected)

    real(wp),allocatable :: Q(:,:)        !> (n3, max_cols) orthonormal basis
    real(wp),allocatable :: gains(:)      !> (max_cols) information gain per step
    real(wp),allocatable :: step_norms(:) !> (max_cols) original step norms

    real(wp) :: gain_threshold = 0.05_wp     !> below this → no new info
    real(wp) :: noise_threshold = 1.0e-7_wp  !> absolute step norm below which
    !> we consider the step noise
    logical  :: frozen = .false.             !> set when freeze is recommended
    integer  :: freeze_count = 0             !> consecutive low-gain steps
    integer  :: freeze_patience = 3          !> freeze after this many consecutive
  end type step_monitor

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

  subroutine monitor_init(mon,nat,max_steps,gain_threshold)
    !*******************************************************************
    !* Initialise the step_monitor.
    !*
    !* mon           : step_monitor object to initialise
    !* nat           : number of atoms
    !* max_steps     : allocation capacity; defaults to 200 if absent
    !* gain_threshold: minimum fractional gain to accept a step
    !*******************************************************************
    type(step_monitor),intent(inout) :: mon
    integer,intent(in)    :: nat
    integer,optional,intent(in)    :: max_steps
    real(wp),optional,intent(in)    :: gain_threshold
    integer :: ms

    ms = 200
    if (present(max_steps)) ms = max_steps

    mon%n3 = 3*nat
    mon%max_cols = ms
    mon%rank = 0
    mon%n_total = 0
    mon%frozen = .false.
    mon%freeze_count = 0

    if (present(gain_threshold)) mon%gain_threshold = gain_threshold

    if (allocated(mon%Q)) deallocate (mon%Q)
    if (allocated(mon%gains)) deallocate (mon%gains)
    if (allocated(mon%step_norms)) deallocate (mon%step_norms)
    allocate (mon%Q(mon%n3,ms),source=0.0_wp)
    allocate (mon%gains(ms),source=0.0_wp)
    allocate (mon%step_norms(ms),source=0.0_wp)
  end subroutine monitor_init

  subroutine monitor_add_step(mon,s_vec,gain,accepted)
    !*******************************************************************
    !* Process a new optimisation step and update the basis Q.
    !*
    !* At each optimisation step, we orthogonalise the new step vector sₖ
    !* against all previous (orthogonalised) steps.  The residual norm
    !* after orthogonalisation is the "information gain" of this step:
    !*
    !*   gain = ‖sₖ − Σᵢ (qᵢᵀsₖ) qᵢ‖  /  ‖sₖ‖
    !*
    !* When gain drops below a threshold, the step adds no new subspace
    !* information → the Hessian should be frozen (further updates corrupt
    !* well-characterised modes with noisy data).
    !* This is equivalent to tracking the rank of the step matrix S via
    !* its thin QR/Gram-Schmidt, but done incrementally.
    !*
    !* mon     : step_monitor object (updated in place)
    !* s_vec   : step vector  Δx = x_{k+1} − x_k  (length 3*nat)
    !* gain    : fractional information gain of this step, ∈ [0, 1]
    !* accepted: .true. if the step expanded the subspace
    !*******************************************************************
    type(step_monitor),intent(inout) :: mon
    real(wp),intent(in)    :: s_vec(:)  ! (n3)
    real(wp),intent(out)   :: gain
    logical,intent(out)   :: accepted

    real(wp),allocatable :: q_new(:)
    real(wp) :: snorm,residual_norm,proj
    integer  :: i

    accepted = .false.
    gain = 0.0_wp
    mon%n_total = mon%n_total+1

    snorm = sqrt(sum(s_vec**2))

    ! ── reject trivially small steps (numerical noise) ───────────────
    if (snorm < mon%noise_threshold) then
      mon%freeze_count = mon%freeze_count+1
      mon%step_norms(mon%n_total) = snorm
      mon%gains(mon%n_total) = 0.0_wp
      return
    end if

    ! ── orthogonalise against existing basis ─────────────────────────
    allocate (q_new(mon%n3))
    q_new = s_vec
    do i = 1,mon%rank
      proj = dot_product(mon%Q(:,i),q_new)
      q_new = q_new-proj*mon%Q(:,i)
    end do

    residual_norm = sqrt(sum(q_new**2))
    gain = residual_norm/snorm

    mon%step_norms(mon%n_total) = snorm
    mon%gains(mon%n_total) = gain

    ! ── accept or reject based on information gain ───────────────────
    if (gain > mon%gain_threshold.and.mon%rank < mon%max_cols) then
      mon%rank = mon%rank+1
      mon%Q(:,mon%rank) = q_new/residual_norm
      accepted = .true.
      mon%freeze_count = 0
    else
      mon%freeze_count = mon%freeze_count+1
    end if

    deallocate (q_new)
  end subroutine monitor_add_step

  logical function monitor_should_freeze(mon)
    !*******************************************************************
    !* Return .true. when freeze_patience consecutive steps have had
    !* low information gain, indicating further Hessian updates are
    !* likely to corrupt rather than improve the approximation.
    !*
    !* mon : step_monitor object
    !*******************************************************************
    type(step_monitor),intent(inout) :: mon

    if (mon%freeze_count >= mon%freeze_patience) mon%frozen = .true.
    monitor_should_freeze = mon%frozen
  end function monitor_should_freeze

  subroutine mode_quality_analysis(modes,n_modes,mon,quality)
    !*******************************************************************
    !* Assess how well each Hessian eigenvector was sampled during opt.
    !*
    !* The quality of mode vᵢ is its squared projection onto the step
    !* subspace:
    !*   quality(i) = ‖Q Qᵀ vᵢ‖²  ∈ [0, 1]
    !*   1.0 ≙ mode thoroughly explored; eigenvalue should be accurate
    !*   0.0 ≙ mode never sampled; eigenvalue is just the initial guess
    !*
    !* modes   : (n3, n_modes) Hessian eigenvectors
    !* n_modes : number of modes
    !* mon     : step_monitor with the orthonormal step subspace Q
    !* quality : (n_modes) sampling quality for each mode, output
    !*******************************************************************
    real(wp),intent(in)  :: modes(:,:)   ! (n3, n_modes)
    integer,intent(in)  :: n_modes
    type(step_monitor),intent(in)  :: mon
    real(wp),intent(out) :: quality(:)   ! (n_modes)

    real(wp),allocatable :: proj_vec(:)
    real(wp) :: proj
    integer  :: i,j

    allocate (proj_vec(mon%n3))

    do i = 1,n_modes
      ! ── compute projection Q Qᵀ vᵢ ──────────────────────────────────
      proj_vec = 0.0_wp
      do j = 1,mon%rank
        proj = dot_product(mon%Q(:,j),modes(:,i))
        proj_vec = proj_vec+proj*mon%Q(:,j)
      end do

      ! ── quality = |Q Qᵀ vᵢ|² (vᵢ is unit-normalised) ───────────────
      quality(i) = sum(proj_vec**2)
      quality(i) = max(0.0_wp,min(1.0_wp,quality(i)))   ! clamp [0,1]
    end do

    deallocate (proj_vec)
  end subroutine mode_quality_analysis

  subroutine selective_hessian_repair(mol,calc,modes,eigenvalues, &
                                      quality,n_modes,q_threshold, &
                                      delta,n_repaired)
    !*******************************************************************
    !* Recompute Hessian eigenvalues for poorly-sampled modes via
    !* central-difference numerical Hessian columns.
    !*
    !* For each mode vᵢ with quality below q_threshold:
    !*   H·vᵢ ≈ [ g(x + δvᵢ) − g(x − δvᵢ) ] / (2δ)
    !*   λᵢ_corrected = vᵢᵀ · (H·vᵢ)
    !* Cost: 2 engrad calls per repaired mode.
    !*
    !* mol         : equilibrium geometry (coord, xyz in Bohr)
    !* calc        : calculation settings (calcdata)
    !* modes       : (n3, n_modes) Hessian eigenvectors
    !* eigenvalues : (n_modes) approx. eigenvalues, corrected in place
    !* quality     : (n_modes) mode quality from mode_quality_analysis
    !* n_modes     : number of modes
    !* q_threshold : quality below this triggers recomputation
    !* delta       : finite-difference displacement in Bohr (e.g. 5e-3)
    !* n_repaired  : number of modes that were recomputed, output
    !*******************************************************************
    type(coord),intent(in)    :: mol
    type(calcdata),intent(inout) :: calc
    real(wp),intent(in)    :: modes(:,:)      ! (n3, n_modes)
    real(wp),intent(inout) :: eigenvalues(:)  ! (n_modes)
    real(wp),intent(in)    :: quality(:)      ! (n_modes)
    integer,intent(in)     :: n_modes
    real(wp),intent(in)    :: q_threshold
    real(wp),intent(in)    :: delta
    integer,intent(out)    :: n_repaired

    type(coord) :: mol_plus,mol_minus
    real(wp),allocatable :: grad_plus(:,:),grad_minus(:,:)  ! (3, nat)
    real(wp),allocatable :: Hv(:),mode_3d(:,:)
    real(wp) :: e_plus,e_minus,lambda_old,lambda_new
    integer  :: i,stat,nat,n3

    nat = mol%nat
    n3 = 3*nat
    n_repaired = 0

    ! ── initialise displaced geometry objects ────────────────────────
    call mol_plus%copy(mol)
    call mol_minus%copy(mol)

    allocate (grad_plus(3,nat),grad_minus(3,nat))
    allocate (Hv(n3),mode_3d(3,nat))

    write (stdout,'(/,a)') ' --- Selective Hessian repair ---'
    write (stdout,'(a,f5.2)') '  Quality threshold: ',q_threshold
    write (stdout,'(a,es9.2,a)') '  FD displacement:   ',delta,' Bohr'

    do i = 1,n_modes

      if (quality(i) >= q_threshold) cycle   ! mode is fine, skip

      n_repaired = n_repaired+1
      lambda_old = eigenvalues(i)

      ! ── reshape mode to (3, nat) for coordinate displacement ─────────
      mode_3d = reshape(modes(:,i), [3,nat])

      ! ── forward displacement ─────────────────────────────────────────
      mol_plus%xyz = mol%xyz+delta*mode_3d
      call engrad(mol_plus,calc,e_plus,grad_plus,stat)

      ! ── backward displacement ────────────────────────────────────────
      mol_minus%xyz = mol%xyz-delta*mode_3d
      call engrad(mol_minus,calc,e_minus,grad_minus,stat)

      ! ── central-difference Hessian–vector product ────────────────────
      Hv = reshape(grad_plus-grad_minus, [n3])/(2.0_wp*delta)

      ! ── corrected eigenvalue:  λ = vᵀ H v ───────────────────────────
      lambda_new = dot_product(modes(:,i),Hv)
      eigenvalues(i) = lambda_new

      write (stdout,'(a,i4,a,f6.3,a,f12.6,a,f12.6,a,f10.6)') &
        '  Mode ',i, &
        '  qual=',quality(i), &
        '  λ_old=',lambda_old, &
        '  λ_new=',lambda_new, &
        '  Δλ=',abs(lambda_new-lambda_old)

    end do

    write (stdout,'(a,i0,a,i0,a)') '  Repaired ',n_repaired,' of ',n_modes,' modes'
    write (stdout,'(a,i0,a)') '  Cost: ',2*n_repaired,' engrad calls'

    call mol_plus%deallocate()
    call mol_minus%deallocate()
    deallocate (grad_plus,grad_minus,Hv,mode_3d)
  end subroutine selective_hessian_repair

! ══════════════════════════════════════════════════════════════════════════════
end module hessian_quality
