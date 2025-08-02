project('edmd', 'c',
  version: '1.0.0',
  default_options: [
    'warning_level=3',
    'optimization=3',
    'c_std=c11'
  ],
  meson_version: '>= 0.56.0'
)

# Add the PseudoRNG subdirectory
subdir('PseudoRNG')

# Math library dependency
math_dep = meson.get_compiler('c').find_library('m', required: true)

# Main executable
edmd_sources = files('mdSingle.c')

edmd_exe = executable(
  'edmd',
  edmd_sources,
  dependencies: [pseudorng_dep, math_dep],
  install: true
)
