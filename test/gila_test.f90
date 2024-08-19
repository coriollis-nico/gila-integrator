program io_test
  !! Test for [[lib_io]] module
  use iso_fortran_env,  &
    only: sp => real32, &
          dp => real64, &
          qp => real128,&
          io => output_unit
  use gila
  implicit none

  type(gila_conditions) :: default_conditions

  character(*), parameter :: separator = &
    &"-----------------------------------------&
    &-----------------------------------------"

! ------------------------------------------------------------------------------------ !

! save_gila_genconditions

print *, "call save_gila_genconditions(gila_grdefault, [output_unit])"
          call save_gila_genconditions(gila_grdefault, io)

print *, "call save_gila_genconditions(default_conditions, [output_unit])"
          call save_gila_genconditions(default_conditions, io)


print *, separator
! gila_friedmann

print *, "gila_friedmann(0.0_qp, 0.0_qp, gila_grdefault)"
print *,  gila_friedmann(0.0_qp, 0.0_qp, gila_grdefault)

print *, "gila_friedmann(0.0_qp, 0.0_qp)"
print *,  gila_friedmann(0.0_qp, 0.0_qp)

end program io_test
