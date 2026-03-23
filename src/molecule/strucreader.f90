!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020-2023 Philipp Pracht
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

!=========================================================================================!
! STRUCRD is a module for reading and writing molecular structures.
!
! The source is organized as follows:
!   0. Variable declarations
!   1. Routines for reading and writing ensemble files/trajectories in the XYZ format
!   2. Routines for reading single structures in various formats
!   3. Routines for writing structures in various formats
!   4. Utility routines mainly used only within the module
!
! Currently supported formats:
!   .xyz (Xmol) files and trajectories (read and write)
!   coord (turbomole) files (read and write)
!   .sdf/.mol files (V2000, read only)
!   .pdb files (in development)
!
!=========================================================================================!
module strucrd
  use molecule_type
  use molecule_type_components
  use molecule_type_ensemble
  use molecule_io
  implicit none
  private
  !> RE-EXPORTS FROM THE ABOVE MODULES
! ══════════════════════════════════════════════════════════════════════════════
  public :: coord     !> coord type
  public :: ensemble  !> ensemble type (sparsely used)
  public :: mollist   !> list of coord objects

!=========================================================================================!
  public :: i2e          !> function to convert atomic number to element symbol
  public :: asym         !> "
  public :: e2i          !> function to convert element symbol into atomic number

  public :: grepenergy
  public :: checkcoordtype
  public :: rdnat       !-- procedure to read number of atoms Nat
  public :: rdcoord     !-- read an input file, determine format automatically
  public :: rdxmol      !-- read a file in the Xmol (.xyz) format specifically
  public :: rdxmolselec !-- read only a certain structure in Xmol file

  public :: wrc0
  public :: wrcoord
  public :: wrxyz
  public :: wrsdf

  public :: xyz2coord
  public :: coord2xyz
  public :: rdensembleparam   !-- read Nat and Nall for a XYZ trajectory
  public :: rdensemble        !-- read a XYZ trajectory
  public :: wrensemble

  public :: coordline
  public :: get_atlist
  public :: sumform

!=========================================================================================!
!=========================================================================================!
contains  !> MODULE PROCEDURES START HERE

end module strucrd
