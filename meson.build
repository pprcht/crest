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

project(
  'crest',
  ['c', 'fortran'],
  version         : '3.0.3',
  license         : 'LGPL-3.0-or-later',
  meson_version   : '>=0.57.0',
  default_options : [
    'buildtype=release',
    'c_std=gnu11',
    'warning_level=0',
  ],
)

# ── Compiler flags, static link args, metadata template, inc_dirs ─────────────
# Sets: fc, cc, fc_id, cc_id, fcid_str, ccid_str, static_build, metadata_conf
subdir('config')

# inc_dirs is defined here rather than in config/meson.build because
# include_directories() rejects absolute paths inside the source/build tree,
# and meson.current_build_dir() only resolves to the *root* build dir when
# called from the root meson.build (subdir() shifts the current dir).
# '.' pairs the source root with its build counterpart automatically,
# which is where configure_file() writes crest_metadata.fh.
inc_dirs = include_directories('include', '.')

# ═══════════════════════════════════════════════════════════════════════════════
#  OpenMP
# ═══════════════════════════════════════════════════════════════════════════════
# ifx (intel-llvm) triggers a Meson bug in dependency('openmp') because its
# preprocessor output lacks the delimiter strings Meson searches for.
# Work around it by constructing the dependency manually for that compiler.
if fc_id == 'intel-llvm'
  if get_option('openmp')
    omp_dep = declare_dependency(
      compile_args : ['-qopenmp'],
      link_args    : ['-qopenmp'],
    )
  else
    omp_dep = declare_dependency()
  endif
else
  omp_dep = dependency('openmp',
    required : get_option('openmp'),
    language : 'fortran',
  )
endif
if omp_dep.found()
  add_project_arguments('-DWITH_OMP', language : ['c', 'fortran'])
endif

_omp_link_dep = (omp_dep.found() and not static_build) ? [omp_dep] : []

# ═══════════════════════════════════════════════════════════════════════════════
#  LAPACK / BLAS
# ═══════════════════════════════════════════════════════════════════════════════
#
# Cross-compiler compatibility matrix
# ─────────────────────────────────────────────────────────────────────────────
#  Fortran compiler  │  Recommended provider  │  MKL threading layer
#  ──────────────────┼────────────────────────┼──────────────────────
#  gfortran          │  OpenBLAS (default)    │  mkl_gnu_thread
#  ifort / ifx       │  MKL  ← natural match  │  mkl_intel_thread
#  ifort/ifx + gcc C │  MKL                   │  mkl_intel_thread (Fortran wins)
# ─────────────────────────────────────────────────────────────────────────────
# The Fortran compiler drives the final link step and therefore owns the OpenMP
# runtime.  Choose the MKL threading layer to match the *Fortran* compiler.

lapack_opt    = get_option('lapack')
with_mkl      = false
with_openblas = false
lapack_dep    = []
blas_dep      = []

_prefer_mkl = (fc_id in ['intel', 'intel-llvm']) and (lapack_opt == 'auto')

# ── MKL ───────────────────────────────────────────────────────────────────────
if lapack_opt == 'mkl' or _prefer_mkl
  _mkl = dependency('mkl-sdl', required : false, static : static_build)

  if not _mkl.found()
    _thread_lib = (fc_id == 'gcc') ? 'mkl_gnu_thread' : 'mkl_intel_thread'
    _omp_lib    = (fc_id == 'gcc') ? 'gomp'           : 'iomp5'
    _mkl_parts  = []
    foreach lib : ['mkl_intel_lp64', _thread_lib, 'mkl_core', _omp_lib, 'm', 'dl']
      _l = fc.find_library(lib, required : false, static : static_build)
      if _l.found()
        _mkl_parts += [_l]
      endif
    endforeach
    if _mkl_parts.length() > 0
      _mkl = declare_dependency(dependencies : _mkl_parts)
    endif
  endif

  if _mkl.found()
    lapack_dep = _mkl
    blas_dep   = _mkl
    with_mkl   = true
    message('LAPACK/BLAS provider: Intel MKL')
  elif lapack_opt == 'mkl'
    error('MKL requested but not found. Set PKG_CONFIG_PATH or source setvars.sh.')
  else
    lapack_opt = 'auto'
  endif
endif

# ── OpenBLAS ──────────────────────────────────────────────────────────────────
if lapack_opt in ['auto', 'openblas'] and not with_mkl
  _opa = dependency('openblas', required : false, static : static_build)
  if not _opa.found()
    _opa = dependency('blas', required : false, static : static_build)
  endif
  if _opa.found()
    lapack_dep    = _opa
    blas_dep      = _opa
    with_openblas = true
    message('LAPACK/BLAS provider: OpenBLAS')
  elif lapack_opt == 'openblas'
    error('OpenBLAS requested but not found. Install libopenblas-dev.')
  endif
endif

# ── Netlib reference ──────────────────────────────────────────────────────────
if lapack_opt in ['auto', 'netlib'] and not with_mkl and not with_openblas
  _blas   = dependency('blas',   required : false, static : static_build)
  _lapack = dependency('lapack', required : false, static : static_build)
  if _blas.found() and _lapack.found()
    blas_dep   = _blas
    lapack_dep = _lapack
    message('LAPACK/BLAS provider: netlib reference')
  else
    _blas   = fc.find_library('blas',   required : false, static : static_build)
    _lapack = fc.find_library('lapack', required : false, static : static_build)
    if _blas.found() and _lapack.found()
      blas_dep   = declare_dependency(dependencies : [_blas])
      lapack_dep = declare_dependency(dependencies : [_lapack])
      message('LAPACK/BLAS provider: system libraries (direct probe)')
    elif lapack_opt == 'netlib'
      error('Netlib BLAS/LAPACK not found.')
    else
      error(
        'No LAPACK/BLAS found. Install libopenblas-dev, intel-mkl, or '
        + 'liblapack-dev+libblas-dev, or use -Dlapack=custom.',
      )
    endif
  endif
endif

# ── Custom ────────────────────────────────────────────────────────────────────
if lapack_opt == 'custom'
  _ll = get_option('lapack_libs')
  _bl = get_option('blas_libs')
  if _ll.length() == 0 and _bl.length() == 0
    error('lapack=custom requires lapack_libs and/or blas_libs.')
  endif
  _deps = []
  foreach lib : _ll + _bl
    _deps += [fc.find_library(lib, required : true, static : static_build)]
  endforeach
  lapack_dep = declare_dependency(dependencies : _deps)
  blas_dep   = lapack_dep
endif

if with_mkl
  add_project_arguments('-DWITH_MKL',      language : ['c', 'fortran'])
endif
if with_openblas
  add_project_arguments('-DWITH_OPENBLAS', language : ['c', 'fortran'])
endif

# ── Static-build guard: pin BLAS/LAPACK for subprojects ───────────────────────
# When building statically, subprojects must not resolve BLAS/LAPACK to shared
# libraries.  meson.override_dependency() is global: all subsequent dependency()
# calls in this project AND every subproject get the static dep we already found.
if static_build
  meson.override_dependency('lapack', lapack_dep)
  meson.override_dependency('blas',   blas_dep)
  if with_openblas
    meson.override_dependency('openblas', blas_dep)
  endif
endif

# libgfortran.a (Fortran runtime) references quadmath_* symbols; the GFortran
# driver would add -lquadmath automatically, but meson uses the C linker driver
# for the final link step, so we must add it explicitly for static GCC builds.
_quadmath_dep = []
if static_build and fc_id == 'gcc'
  _qm = fc.find_library('quadmath', required : false, static : true)
  if _qm.found()
    _quadmath_dep = [_qm]
  endif
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Optional subproject / external library dependencies
# ═══════════════════════════════════════════════════════════════════════════════
tomlf_dep = dependency('toml-f',
  version  : '>=0.2.4',
  fallback : ['toml-f', 'tomlf_dep'],
  required : get_option('toml-f'),
  static   : static_build,
)
with_tomlf = tomlf_dep.found()
if with_tomlf
  add_project_arguments('-DWITH_TOMLF', language : ['c', 'fortran'])
endif

tblite_dep = dependency('tblite',
  version  : '>=0.3.0',
  fallback : ['tblite', 'tblite_dep'],
  required : get_option('tblite'),
  static   : static_build,
  default_options: ['api=false'],
)
with_tblite = tblite_dep.found()
if with_tblite
  add_project_arguments('-DWITH_TBLITE', language : ['c', 'fortran'])
endif

gfnff_dep = dependency('gfnff',
  fallback : ['gfnff', 'gfnff_dep'],
  required : get_option('gfnff'),
  static   : static_build,
  default_options : ['tests=false'],
)
with_gfnff = gfnff_dep.found()
if with_gfnff
  add_project_arguments('-DWITH_GFNFF', language : ['c', 'fortran'])
endif

gfn0_dep = dependency('gfn0',
  fallback : ['gfn0', 'gfn0_dep'],
  required : get_option('gfn0'),
  static   : static_build,
)
with_gfn0 = gfn0_dep.found()
if with_gfn0
  add_project_arguments('-DWITH_GFN0', language : ['c', 'fortran'])
endif

libpvol_dep = dependency('libpvol',
  fallback : ['pvol', 'libpvol_dep'],
  required : get_option('libpvol'),
  static   : static_build,
  default_options : ['tests=false'],
)
with_libpvol = libpvol_dep.found()
if with_libpvol
  add_project_arguments('-DWITH_LIBPVOL', language : ['c', 'fortran'])
endif

lwoniom_dep = dependency('lwoniom',
  fallback : ['lwoniom', 'lwoniom_dep'],
  required : get_option('lwoniom'),
  static   : static_build,
)
with_lwoniom = lwoniom_dep.found()
if with_lwoniom
  add_project_arguments('-DWITH_LWONIOM', language : ['c', 'fortran'])
endif

fmlip_dep = dependency('fmlip_relay',
  fallback : ['fmlip_relay', 'fmlip_relay_dep'],
  required : get_option('fmlip-relay'),
  static   : static_build,
)
with_fmlip = fmlip_dep.found()
if with_fmlip
  add_project_arguments('-DWITH_FMLIP_RELAY', language : ['c', 'fortran'])
endif

if get_option('tests')
  testdrive_dep = dependency('test-drive',
    version  : '>=0.4.0',
    fallback : ['test-drive', 'testdrive_dep'],
    required : true,
    static   : static_build,
  )
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Metadata header — finalise and generate
# ═══════════════════════════════════════════════════════════════════════════════
# config/meson.build populated metadata_conf with system/compiler info but
# deferred the with_* feature flags (not yet known at that point).
# Patch them in now that all deps have been resolved.
metadata_conf.set('tomlfvar',   with_tomlf.to_string())
metadata_conf.set('gfn0var',    with_gfn0.to_string())
metadata_conf.set('gfnffvar',   with_gfnff.to_string())
metadata_conf.set('tblitevar',  with_tblite.to_string())
metadata_conf.set('libpvolvar', with_libpvol.to_string())
metadata_conf.set('lwoniomvar', with_lwoniom.to_string())

configure_file(
  input         : 'assets/template/metadata.f90',
  output        : 'crest_metadata.fh',
  configuration : metadata_conf,
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Source collection
# ═══════════════════════════════════════════════════════════════════════════════
srcs = []
prog = []
subdir('src')

# ═══════════════════════════════════════════════════════════════════════════════
#  Build targets
# ═══════════════════════════════════════════════════════════════════════════════
_optional_deps = []
foreach d : [tblite_dep, gfn0_dep, gfnff_dep, libpvol_dep, lwoniom_dep, tomlf_dep, fmlip_dep]
  if d.found()
    _optional_deps += [d]
  endif
endforeach

lib_crest = static_library(
  'crest',
  sources             : srcs,
  include_directories : inc_dirs,
  dependencies        : [lapack_dep, blas_dep] + _omp_link_dep + _optional_deps + _quadmath_dep,
  install             : false,
  pic                 : true,
)

executable(
  'crest',
  sources             : prog,
  include_directories : inc_dirs,
  link_with           : lib_crest,
  dependencies        : [lapack_dep, blas_dep] + _omp_link_dep + _optional_deps + _quadmath_dep,
  link_language       : 'fortran',
  install             : true,
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Docs, tests, install
# ═══════════════════════════════════════════════════════════════════════════════
subdir('docs')

if get_option('tests')
  subdir('test')
endif

install_data(
  'COPYING', 'COPYING.LESSER',
  install_dir : get_option('datadir') / 'licenses' / meson.project_name(),
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════════════════════
summary({
  'Fortran compiler' : '@0@ @1@'.format(fcid_str, fc.version()),
  'C compiler'       : '@0@ @1@'.format(ccid_str, cc.version()),
  'Build type'       : get_option('buildtype'),
  'Static binary'    : static_build,
  'OpenMP'           : omp_dep.found(),
  'MKL'              : with_mkl,
  'OpenBLAS'         : with_openblas,
  'tblite'           : with_tblite,
  'toml-f'           : with_tomlf,
  'GFN-FF'           : with_gfnff,
  'GFN0-xTB'         : with_gfn0,
  'libpvol'          : with_libpvol,
  'lwONIOM'          : with_lwoniom,
  'fmlip-relay'      : with_fmlip,
  'Unit tests'       : get_option('tests'),
}, section : 'crest build configuration')
