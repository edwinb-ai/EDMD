# PseudoRNG library
pseudorng_sources = files(
  'xoshiro256plus.c',
  'splitmix64.c'
)

pseudorng_lib = static_library(
  'PseudoRNG',
  pseudorng_sources,
  install: false
)

pseudorng_dep = declare_dependency(
  link_with: pseudorng_lib,
  include_directories: include_directories('.')
)
