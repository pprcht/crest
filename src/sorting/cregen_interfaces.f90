
!=========================================================================================!  
!=========================================================================================!  
!> Interfaces for use CREGEN (and related)
!=========================================================================================!  
!=========================================================================================!  

module cregen_interface
!*******************************************************
!* module to load an interface to the newcregen routine
!* mandatory to handle the optional input arguments
!*******************************************************
  use unionize_module
  implicit none
  interface
    subroutine newcregen(env,quickset,infile)
      use crest_parameters
      use crest_data
      use crest_restartlog
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      integer,intent(in),optional :: quickset
      character(len=*),intent(in),optional :: infile
    end subroutine newcregen

    subroutine cregen_irmsd_all(nall,structures,printlvl,iinversion)
      use strucrd
      implicit none
      !> INPUT
      integer,intent(in) :: nall
      type(coord),intent(inout),target :: structures(nall)
      integer,intent(in),optional :: printlvl
      integer,intent(in),optional :: iinversion
    end subroutine cregen_irmsd_all

    subroutine cregen_irmsd_sort(env,nall,structures,groups,allcanon,printlvl)
      use crest_data
      use strucrd
      implicit none
      !> INPUT
      type(systemdata),intent(inout) :: env
      integer,intent(in) :: nall
      type(coord),intent(inout),target :: structures(nall)
      integer,intent(inout) :: groups(nall)
      logical,intent(in),optional :: allcanon
      integer,intent(in),optional :: printlvl
    end subroutine cregen_irmsd_sort

  end interface
!>--- Additional Related RE-EXPORTS
  public :: unionizeEnsembles
end module cregen_interface

!=========================================================================================! 
!=========================================================================================! 
!> Interfaces for routines used WITHIN CREGEN
!=========================================================================================! 
!=========================================================================================! 

module cregen_subroutines
!*************************************
!* interfaces for cregen subroutines
!*************************************
  implicit none
  interface
    subroutine discardbroken(ch,env,topocheck,structures,newnall)
      use crest_data
      use strucrd
      use cregen_utils
      type(systemdata),intent(in) :: env
      integer,intent(in) :: ch
      logical,intent(in) :: topocheck
      type(coord),intent(inout),allocatable,target :: structures(:)
      integer,intent(out) :: newnall
    end subroutine discardbroken

  end interface
end module cregen_subroutines
