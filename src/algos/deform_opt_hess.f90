subroutine deform_opt_hess(calc,mol)  
    use crest_calculator
    use strucrd
    use irmsd_module
    use bh_step_module
    use crest_parameters
    use optimize_module
    implicit none
    type(calcdata),intent(inout) :: calc
    type(coord),intent(in) :: mol
    type(coord) :: molnew,mol_reopt
    real(wp) :: energy,stepsize,rmsdval
    real(wp) :: grad(3,mol%nat) 
    logical :: pr,wr 
    integer :: io

    if (allocated(calc%chess)) deallocate(calc%chess)

    !allocate (calc%chess)
    !call calc%chess%alloc(mol%nat,calc%hu_steps,calc%hguess,calc%initialize_hr_type, calc%hr_hu_type) !Maybe in future just reset the cash here, or if this works we only call this here

    molnew=mol

    stepsize = 0.25_wp

    call take_fixed_stepsize_cart(molnew,stepsize,calc)

    pr = .true.
    wr = .true.

    call optimize_geometry(molnew,mol_reopt,calc,energy,grad,pr,wr,io)

    rmsdval = rmsd(mol,mol_reopt)

    write(stdout,*) 'VALUE OF RMSD FOR REOPTIMISED STRUCTURE',rmsdval


end subroutine deform_opt_hess