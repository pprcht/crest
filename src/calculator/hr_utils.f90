module hr_utils
  use iso_fortran_env,only:wp => real64
  use crest_calculator
  use crest_parameters
  implicit none
  private

  public hr_initialize_hessian

contains

subroutine hr_initialize_hessian(calc,at)
  type(calcdata),intent(inout) :: calc
  type(calcdata) :: newcalc
  type(calculation_settings) :: clevel
  integer :: k,idx,io, nat3
  integer, intent(in) :: at(:)
  
  nat3 = 3*calc%chess%natm
  idx = minloc(calc%chess%order,1) !>gives location of first geometry that is saved

    !>initialize_type: 0 for scaled identity, 1 for gfnff guess, 2 for gfn2 guess

select case (calc%chess%initialize_type)
case(0)
    calc%chess%hguess_mat = 0.0_wp
    do k = 1,nat3
        calc%chess%hguess_mat(k,k) = calc%chess%hguess
    end do
case(1)
    call clevel%create('gfnff', chrg=calc%calcs(1)%chrg, uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
    call newcalc%add(clevel)
    call numhess1(calc%chess%natm,at,calc%chess%coords(idx,:,:),newcalc,calc%chess%hguess_mat(:,:),io)      
case(2)
    call clevel%create('gfn2', chrg=calc%calcs(1)%chrg, uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
    call newcalc%add(clevel)
    call numhess1(calc%chess%natm,at,calc%chess%coords(idx,:,:),newcalc,calc%chess%hguess_mat(:,:),io)
end select
end subroutine hr_initialize_hessian

end module hr_utils