# This file is part of crest.
# SPDX-Identifier: LGPL-3.0-or-later
#
# crest is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# crest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with crest.  If not, see <https://www.gnu.org/licenses/>.

set(_lib "lammps-mace")
set(_pkg "LAMMPS-MACE")
set(_url "https://github.com/ACEsuit/lammps")


#########################################################################################
set(package "lammps")
string(TOLOWER "${package}" _pkg_lc)
string(TOUPPER "${package}" _pkg_uc)

if(NOT DEFINED "${_pkg_uc}_SUBPROJECT")
  set("_${_pkg_uc}_SUBPROJECT")
  set("${_pkg_uc}_SUBPROJECT" "subprojects/${package}/cmake")
endif()
set("${_pkg_uc}_SOURCE_DIR" "${PROJECT_SOURCE_DIR}/${${_pkg_uc}_SUBPROJECT}")
set("${_pkg_uc}_BINARY_DIR" "${PROJECT_BINARY_DIR}/${${_pkg_uc}_SUBPROJECT}")
if(EXISTS "${${_pkg_uc}_SOURCE_DIR}/CMakeLists.txt")
  message(STATUS "Include ${package} from ${${_pkg_uc}_SUBPROJECT}")
#  add_subdirectory(
#    "${${_pkg_uc}_SOURCE_DIR}"
#    "${${_pkg_uc}_BINARY_DIR}"
#  )


  # MACE requires libtorch
  # The FetchContent can be tricky, I needed to upgrade my CMake to 3.28
  # and the following policy needs to be set
  if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
		cmake_policy(SET CMP0135 NEW)
	endif()
  include(FetchContent)
  FetchContent_Declare(libtorch URL "https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.13.0%2Bcpu.zip")
  FetchContent_MakeAvailable(libtorch)
  set(CMAKE_PREFIX_PATH "${libtorch_SOURCE_DIR};${CMAKE_PREFIX_PATH}")
  message(STATUS "${libtorch_SOURCE_DIR}")

  # LAMMPS build options
  OPTION(BUILD_MPI "MPI LAMMPS build" OFF)
  #OPTION(PKG_MOLECULE "" ON)
  OPTION(BUILD_OMP "" OFF)
  OPTION(PKG_OPENMP "" OFF)
  OPTION(PKG_ML-MACE "MACE option" ON)
  
  add_library("${package}::${package}" INTERFACE IMPORTED)
  target_link_libraries("${package}::${package}" INTERFACE "${package}")

  # We need the module directory in the subproject before we finish the configure stage
  if(NOT EXISTS "${${_pkg_uc}_BINARY_DIR}/include")
    make_directory("${${_pkg_uc}_BINARY_DIR}/include")
  endif()

else()
message(FATAL_ERROR "LAMMPS + MACE could not be found! Please place the source in subprojects/lammps-mace")
endif()

#########################################################################################

set(found FALSE)
if(TARGET "${package}::${package}")
  set (found TRUE)
endif()
message(STATUS "Found LAMMPS + MACE : ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
unset(package)
