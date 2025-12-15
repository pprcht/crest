# This file is part of crest.
# SPDX-Identifier: LGPL-3.0-or-later

project(
  'crest',
  ['fortran', 'c'],
  version: '3.0.3',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.63',
  default_options: [
    'buildtype=debugoptimized',
    'warning_level=0',
    # IMPORTANT:
    # Do NOT force -static here. Linking mode should be controlled via options
    # (e.g. -Dstatic_exe / -Dlink_mode) implemented in config/meson.build.
    #
    # Also do NOT force default_library here. Users can choose via:
    #   -Ddefault_library=static|shared
  ],
)

# Install only when this is the main project.
# Subprojects should generally not install by default.
do_install = not meson.is_subproject()

# =========================
# Configuration (toolchain, deps, options)
# =========================
exe_deps = []
config = configuration_data()

# config/meson.build is responsible for:
# - enforcing compiler consistency (FC/CC)
# - adding compiler/link flags
# - selecting OpenMP
# - selecting BLAS/LAPACK backend
# - populating `exe_deps`
# - populating `config` (configuration_data)
subdir('config')

# Generate configured metadata include
configure_file(
  input: files('assets/template/metadata.f90'),
  output: 'crest_metadata.fh',
  configuration: config,
)

# =========================
# Sources
# =========================
prog = []
srcs = []
subdir('src')

# =========================
# Targets
# =========================
crest_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: exe_deps,
  install: do_install,  
)

crest_inc = include_directories('.')  # safest baseline; add 'include' if you have one

crest_dep = declare_dependency(
  link_with: crest_lib,
  include_directories: crest_inc,
  dependencies: exe_deps,
)

# Allow `dependency('crest')` when used as a subproject
meson.override_dependency(meson.project_name(), crest_dep)

crest_exe = executable(
  meson.project_name(),
  sources: prog,
  dependencies: crest_dep,
  install: do_install,
  link_language: 'fortran',
)

summary(
  {
    'Install': do_install,
    'default_library': get_option('default_library'),
    'OpenMP option': get_option('openmp'),
    'LA backend': get_option('la_backend'),
  },
  section: 'crest',
)

