project('mdgrow', 'c',
  version: '1.0.0',
  default_options: [
    'warning_level=3',
    'optimization=3',
    'c_std=c11'
  ],
  meson_version: '>= 0.56.0'
)

# Math library dependency
math_dep = meson.get_compiler('c').find_library('m', required: true)

# Source files - only mdgrow.c since it includes mt19937ar.c directly
mdgrow_sources = files('mdgrow.c')

# Main executable
mdgrow_exe = executable(
  'mdgrow',
  mdgrow_sources,
  dependencies: [math_dep],
  c_args: ['-march=native', '-funroll-loops'],
  install: true
)
