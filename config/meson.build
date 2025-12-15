# This file is part of crest.
# SPDX-Identifier: LGPL-3.0-or-later

cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')

cc_id = cc.get_id()
fc_id = fc.get_id()

# ------------------------------------------------------------ 
# Enforce toolchain consistency (no mixing)                    
# ------------------------------------------------------------ 
# Normalize compiler ids into families
cc_family = cc_id
if cc_id.startswith('intel')
  cc_family = 'intel'
elif cc_id == 'gcc'
  cc_family = 'gcc'
elif cc_id == 'clang'
  cc_family = 'clang'
endif

fc_family = fc_id
if fc_id.startswith('intel')
  fc_family = 'intel'
elif fc_id == 'gcc'
  fc_family = 'gcc'
elif fc_id == 'clang'
  fc_family = 'clang'
endif

if cc_family != fc_family
  error('Mixed toolchains are not supported: CC=' + cc_id + ' (family ' + cc_family +
        ') but FC=' + fc_id + ' (family ' + fc_family + ').')
endif

toolchain = fc_family
message('Toolchain family: ' + toolchain + ' (CC=' + cc_id + ', FC=' + fc_id + ')')

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------
openmp_opt   = get_option('openmp')     # 'auto'|'true'|'false'
link_mode    = get_option('link_mode')  # 'default'|'prefer-static'|'fully-static'
mkl_link     = get_option('mkl_link')   # 'auto'|'rt'|'static'
custom_libs  = get_option('custom_libraries')
la_backend = get_option('la_backend')
if la_backend == 'auto'
  if toolchain == 'gcc'
    la_backend = 'openblas'
  elif toolchain == 'intel' or toolchain == 'intel-cl'
    la_backend = 'mkl'
  else
    la_backend = 'netlib'
  endif
  message('LA backend: auto → ' + la_backend)
endif

# ------------------------------------------------------------
# Link mode handling
# ------------------------------------------------------------
if link_mode == 'prefer-static'
  # Prefer static libraries when available (Meson built-in option).
  # This avoids forcing a fully static executable (which is much more fragile).
  meson.override_option('prefer_static', true)
  message('Link mode: prefer-static (Meson prefer_static=true)')
elif link_mode == 'fully-static'
  message('Link mode: fully-static (-static will be attempted; this may fail on some systems)')
  # Apply to both C and Fortran link steps (final link is Fortran, but keep consistent).
  if cc.has_link_argument('-static')
    add_project_link_arguments('-static', language: 'c')
  endif
  if fc.has_link_argument('-static')
    add_project_link_arguments('-static', language: 'fortran')
  endif
else
  message('Link mode: default')
endif

# ------------------------------------------------------------
# OpenMP handling (auto/true/false)
# ------------------------------------------------------------
omp_deps = []
openmp_enabled = false

if openmp_opt != 'false'
  omp_required = (openmp_opt == 'true')

  # Language-aware OpenMP deps. For crest (Fortran main), Fortran is the important one.
  omp_f = dependency('openmp', language: 'fortran', required: omp_required)
  if omp_f.found()
    omp_deps += omp_f
    openmp_enabled = true
  endif

  # If you compile any C sources that use OpenMP, this helps too.
  # Not required unless you know you need it.
  omp_c = dependency('openmp', language: 'c', required: false)
  if omp_c.found()
    omp_deps += omp_c
    openmp_enabled = true
  endif

  if openmp_opt == 'auto' and not openmp_enabled
    message('OpenMP: not found (auto → continuing without OpenMP)')
  elif openmp_opt == 'true'
    message('OpenMP: enabled (required)')
  else
    message('OpenMP: enabled')
  endif
else
  message('OpenMP: disabled by option')
endif

omp_dep = declare_dependency(dependencies: omp_deps)

# ------------------------------------------------------------
# Linear algebra backend (BLAS/LAPACK)
# Produces `la_dep` that can be added to dependencies.
# ------------------------------------------------------------
la_dep = declare_dependency()

threads_dep = dependency('threads')

if la_backend == 'none'
  message('LA backend: none')
elif la_backend == 'netlib'
  message('LA backend: netlib (system blas+lapack)')
  blas = dependency('blas', required: true)
  lapack = dependency('lapack', required: true)
  la_dep = declare_dependency(dependencies: [blas, lapack])

elif la_backend == 'openblas'
  message('LA backend: OpenBLAS')
  # Prefer pkg-config/cmake detection if available
  openblas = dependency('openblas', required: false)
  if openblas.found()
    # Many distros provide OpenBLAS with LAPACK symbols, but not all.
    # Try LAPACK as well; if it's already included, this is usually harmless.
    lapack = dependency('lapack', required: false)
    la_dep = declare_dependency(dependencies: lapack.found() ? [openblas, lapack] : [openblas])
  else
    # Fallback to plain library lookup
    lib_ob = cc.find_library('openblas', required: true)
    lapack = dependency('lapack', required: false)
    la_dep = declare_dependency(dependencies: lapack.found() ? [lib_ob, lapack] : [lib_ob])
  endif


elif la_backend == 'custom'
  message('LA backend: custom')
  # Allow passing raw -l... flags or full paths via custom_libraries option
  if custom_libs.length() == 0
    error('la_backend=custom but custom_libraries is empty.')
  endif
  la_dep = declare_dependency(link_args: custom_libs)

elif la_backend == 'mkl'
  message('LA backend: MKL (mkl_link=' + mkl_link + ')')

  # Common system libs MKL often needs
  sys_libs = [threads_dep]
  lib_dl = cc.find_library('dl', required: false)
  if lib_dl.found()
    sys_libs += lib_dl
  endif

  # We reuse the global OpenMP selection:
  # - if openmp_enabled: MKL threaded variant (and OpenMP runtime already handled by omp_dep)
  # - else: MKL sequential
  mkl_is_threaded = openmp_enabled

  # Helper: choose LP64 interface library name based on toolchain.
  # For GNU Fortran, "mkl_gf_lp64" is typical.
  # For Intel, "mkl_intel_lp64" is typical.
  mkl_lp64 = ''
  if toolchain == 'gcc'
    mkl_lp64 = 'mkl_gf_lp64'
  elif toolchain == 'intel' or toolchain == 'intel-cl'
    mkl_lp64 = 'mkl_intel_lp64'
  else
    # If you later want to support clang/nvhpc, extend here.
    warning('MKL selected but toolchain "' + toolchain + '" is not explicitly handled; trying Intel LP64 interface.')
    mkl_lp64 = 'mkl_intel_lp64'
  endif

  if mkl_link == 'rt' or mkl_link == 'auto'
    # Prefer MKL runtime library; simplest and most portable across MKL installs.
    lib_mklrt = cc.find_library('mkl_rt', required: (mkl_link == 'rt'))
    if lib_mklrt.found()
      la_dep = declare_dependency(
        dependencies: sys_libs + [lib_mklrt, omp_dep],
      )
    elif mkl_link == 'rt'
      error('mkl_link=rt requested but libmkl_rt was not found.')
    endif
  endif

  if (mkl_link == 'static') or (mkl_link == 'auto' and not la_dep.found())
    # Static archive linkage (MKL "classic" components).
    # On many systems you need start-group/end-group to resolve circular deps.
    mkl_libs = []
    mkl_libs += cc.find_library(mkl_lp64, required: true)
    mkl_libs += cc.find_library('mkl_core', required: true)

    if mkl_is_threaded
      # Threaded MKL component depends on toolchain
      if toolchain == 'gcc'
        mkl_libs += cc.find_library('mkl_gnu_thread', required: true)
      else
        # Intel toolchain (and default fallback)
        mkl_libs += cc.find_library('mkl_intel_thread', required: true)
      endif
      # OpenMP runtime comes from omp_dep (dependency('openmp'))
      mkl_extra = [omp_dep]
    else
      mkl_libs += cc.find_library('mkl_sequential', required: true)
      mkl_extra = []
    endif

    # Add group flags if supported (mainly relevant for GNU ld)
    grp_args = []
    if fc.has_link_argument('-Wl,--start-group') and fc.has_link_argument('-Wl,--end-group')
      grp_args = ['-Wl,--start-group'] + [] + ['-Wl,--end-group']
    endif

    # Put MKL libs between start/end group if we can
    link_args = []
    if grp_args.length() > 0
      link_args = ['-Wl,--start-group']
      foreach l : mkl_libs
        # Meson "library" objects can be used as dependencies, but not as raw args here.
        # So we keep them as deps and only use group flags as link_args.
      endforeach
      link_args += ['-Wl,--end-group']
    endif

    la_dep = declare_dependency(
      dependencies: sys_libs + mkl_libs + mkl_extra,
      link_args: link_args,
    )
  endif

else
  error('Unknown la_backend value: ' + la_backend)
endif

# ------------------------------------------------------------
# Project warnings / compile args (kept minimal and portable)
# ------------------------------------------------------------
## Add flags only if the compiler supports them.
#foreach flag : ['-Wall', '-Wextra']
#  if cc.has_argument(flag)
#    add_project_arguments(flag, language: 'c')
#  endif
#endforeach
#
## Fortran warnings are much less standardized; keep conservative.
#foreach flag : ['-Wall']
#  if fc.has_argument(flag)
#    add_project_arguments(flag, language: 'fortran')
#  endif
#endforeach

# ------------------------------------------------------------
# Export deps to parent scope
# ------------------------------------------------------------
# Everything the build needs should be captured in dependencies, not global link args.
exe_deps += [omp_dep, la_dep]


# Decide whether fallback subprojects should be built as static libs
subproj_static = (link_mode == 'prefer-static' or link_mode == 'fully-static')
subproj_opts = []
if subproj_static
  subproj_opts = ['default_library=static']
endif

# xtb methods via tblite
if get_option('WITH_TBLITE')
  add_project_arguments('-DWITH_TBLITE', language: 'fortran')

  tblite_dep = dependency(
    'tblite',
    version: '>=0.5',
    fallback: ['tblite', 'tblite_dep'],
    default_options: subproj_opts + ['api=false'],
  )
  exe_deps += tblite_dep
endif

# TOML-F parser
if get_option('WITH_TOMLF')
 add_project_arguments('-DWITH_TOMLF', language: 'fortran')
 tomlf_dep = dependency(
   'toml-f',
   version: '>=0.2.0',
   fallback: ['toml-f', 'tomlf_dep'],
   default_options: subproj_opts + [],
 )
 exe_deps += tomlf_dep
endif

# GFN0-xTB
if get_option('WITH_GFN0')
 add_project_arguments('-DWITH_GFN0', language: 'fortran')
 gfn0_dep = dependency(
   'gfn0',
   fallback: ['gfn0', 'gfn0_dep'],
   default_options: subproj_opts + ['with_gbsa=true'],
 )
 exe_deps += gfn0_dep
endif

# GFN-FF
if get_option('WITH_GFNFF')
 add_project_arguments('-DWITH_GFNFF', language: 'fortran')
 gfnff_dep = dependency(
   'gfnff',
   fallback: ['gfnff', 'gfnff_dep'],
   default_options: subproj_opts + ['with_gbsa=true'],
 )
 exe_deps += gfnff_dep
endif

# LIBPVOL
if get_option('WITH_LIBPVOL')
 add_project_arguments('-DWITH_LIBPVOL', language: 'fortran')
 libpvol_dep = dependency(
   'libpvol',
   fallback: ['libpvol', 'libpvol_dep'],
   default_options: subproj_opts + [],
 )
 exe_deps += libpvol_dep
endif


# lwONIOM
if get_option('WITH_LWONIOM')
 add_project_arguments('-DWITH_LWONIOM', language: 'fortran')
 lwoniom_dep = dependency(
   'lwoniom',
   fallback: ['lwoniom', 'lwoniom_dep'],
   default_options: subproj_opts + [],
 )
 exe_deps += lwoniom_dep
endif



## ========================================= ##
## populate the data for crest_metadata.fh
## ========================================= ##  
commit = 'unknown commit'
git = find_program('git', required: false)
if git.found()
  git_commit = run_command(git, 'show', '-s', '--format=%h',check:true)
  if git_commit.returncode() == 0
    commit = git_commit.stdout().strip()
  endif
endif
# create configuration data
config = configuration_data({
  'name': meson.project_name(),
  'version': meson.project_version(),
  'description': 'Conformer Rotamer Ensemble Sampling Tool',
  'commit': commit,
  'date': run_command('date',check:true).stdout().     strip(),
  'author': run_command('id','-u','-n', check:true).stdout().strip(),
  'origin': run_command('hostname', check:true).stdout().strip(),
  'fcid': fc.get_id(),
  'fcver': fc.version(),
  'ccid': cc.get_id(),
  'ccver': cc.version(),
  'bsystem': 'meson '+meson.version(),
  'tomlfvar': get_option('WITH_TOMLF') ? 'true' : 'false',
  'gfn0var': get_option('WITH_GFN0') ? 'true' : 'false',
  'gfnffvar': get_option('WITH_GFNFF') ? 'true' : 'false', 
  'tblitevar': get_option('WITH_TBLITE') ? 'true' : 'false',
  'libpvolvar': get_option('WITH_LIBPVOL') ? 'true' : 'false',
  'lwoniomvar': get_option('WITH_LWONIOM') ? 'true' : 'false',
})

# Optional: a tiny summary here is helpful while iterating
message('Config summary: la_backend=' + la_backend +
        ', openmp=' + openmp_opt +
        ', openmp_enabled=' + (openmp_enabled ? 'true' : 'false') +
        ', link_mode=' + link_mode)

