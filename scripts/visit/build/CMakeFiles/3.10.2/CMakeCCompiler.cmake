set(CMAKE_C_COMPILER "/opt/cray/pe/craype/2.7.0/bin/cc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "9.3.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "CrayPrgEnv")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "/usr/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "/usr/bin/gcc-ranlib")
set(CMAKE_LINKER "/apps/daint/UES/xalt/xalt2/software/xalt/2.8.10/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "AtpSigHandler;rca;sci_gnu_82_mpi;sci_gnu_82;fftw3_mpi;fftw3_threads;fftw3;fftw3f_mpi;fftw3f_threads;fftw3f;hdf5_hl;hdf5;mpich_gnu_82;gfortran;quadmath;mvec;m;pthread;gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/opt/cray/pe/libsci/20.06.1/GNU/8.1/x86_64/lib;/opt/cray/pe/fftw/3.3.8.7/broadwell/lib;/opt/cray/dmapp/default/lib64;/opt/cray/pe/mpt/7.7.15/gni/mpich-gnu/8.2/lib;/opt/cray/pe/hdf5/1.12.0.0/gnu/8.2/lib;/opt/cray/rca/2.2.20-7.0.2.1_2.27__g8e3fb5b.ari/lib64;/opt/cray/pe/atp/3.7.4/lib;/apps/daint/UES/xalt/xalt2/software/xalt/2.8.10/lib64;/opt/gcc/9.3.0/snos/lib/gcc/x86_64-suse-linux/9.3.0;/opt/gcc/9.3.0/snos/lib64;/lib64;/usr/lib64;/apps/daint/UES/jenkins/7.0.UP02/mc/easybuild/software/GSL/2.5-CrayIntel-20.08/lib;/apps/daint/UES/xalt/xalt2/software/xalt/2.8.10/lib;/opt/gcc/9.3.0/snos/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
