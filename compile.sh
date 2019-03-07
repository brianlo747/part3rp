
gfortran -c -g -O0 microphysics_constants.F90
gfortran -c -g -O0 microphysics_common.F90  
gfortran -c -g -O0 mphys_with_ice.F90

gfortran -o function_tests -g -O0 function_tests.F90 *.o 
