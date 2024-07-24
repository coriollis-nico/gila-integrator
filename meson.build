project(
  'GILA Integrator', 'fortran',
  default_options :
    ['buildtype=debug',
    'warning_level=everything',
    'fortran_std=f2008'],
  license : 'MIT',
  license_files : 'LICENSE.txt'
  )

if meson.get_compiler('fortran').has_argument('-march=native') == true
  add_project_arguments('-march=native', language: 'fortran', native : true)
endif

# lib names
gila_src = 'src/lib/gila.f90'
finite_diff_src = 'src/lib/finite_diff.f90'
lib_io_src = 'src/lib/lib_io.f90'

# test names
gila_test_src = 'src/test/gila_test.f90'
io_test_src = 'src/test/io_test.f90'
fd_test_src = 'src/test/fd_test.f90'

# main names
gr_int_src = 'src/main/gr_integration.f90'
curve_sandwich_src = 'src/main/curve_sandwich.f90'

# Tests
gila_test = executable('gila_test', files(gila_test_src, gila_src, finite_diff_src))
test('gila_test', gila_test)

io_test = executable('io_test', files(io_test_src, lib_io_src))
test('lib_io_test', io_test)

fd_test = executable('fd_test', files(fd_test_src, finite_diff_src))
test('fd_test', fd_test)

# Main
executable('gr_integration', files(gr_int_src, gila_src, lib_io_src, finite_diff_src))
executable('curve_sandwich', files(curve_sandwich_src, gila_src, lib_io_src, finite_diff_src))
