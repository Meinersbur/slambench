cmake_minimum_required(VERSION 2.8)

set(appname kfusion)

find_package(TooN   REQUIRED)
find_package(OpenCL )
find_package(CUDA   )
find_package(GLUT   )
find_package(OpenGL )
find_package(OpenNI2)
find_package(PENCIL )

find_package(PkgConfig)
find_package(Qt5Core)
find_package(Qt5Widgets)
find_package(Qt5Gui)
find_package(Qt5OpenGL)
find_package(Qt5PrintSupport)

include_directories(./include ${TOON_INCLUDE_PATHS} ./thirdparty)


 # ---- PREPARE COMMON DEPENDENCIES  ------------ 
if(APPLE)
    set(common_libraries     stdc++ pthread)
else(APPLE)
    set(common_libraries     rt stdc++ pthread)
endif(APPLE)


set(common_compile_flags "")


 # ----------------- CPP VERSION ----------------- 

message(STATUS "Adding a plain C++ version")
add_library(${appname}-cpp  src/cpp/kernels.cpp)
target_link_libraries(${appname}-cpp   ${common_libraries})	
add_version(${appname} cpp "" "")

 # ----------------- OPENMP VERSION ----------------- 

message(STATUS "Compiling an OpenMP version")
add_library(${appname}-openmp  src/cpp/kernels.cpp)
target_link_libraries(${appname}-openmp   ${common_libraries})	
SET_TARGET_PROPERTIES(${appname}-openmp PROPERTIES COMPILE_FLAGS "-fopenmp")
add_version(${appname} openmp "-fopenmp" "-fopenmp")


 #  ----------------- OCL VERSION ----------------- 
 
if (OPENCL_FOUND) 
  message(STATUS "Adding an OpenCL version")

   include_directories(${OPENCL_INCLUDE_DIRS})
    add_library(${appname}-opencl  src/opencl/kernels.cpp src/opencl/common_opencl.cpp)
    target_link_libraries(${appname}-opencl   ${common_libraries} ${OPENCL_LIBRARIES})	
    add_version(${appname} opencl "" "")
endif(OPENCL_FOUND)


 #  ----------------- CUDA VERSION ----------------- 
if (CUDA_FOUND) 
  message(STATUS "Adding a CUDA version")

set(kfusion_cuda_srcs
	src/cuda/kernels.cu
	thirdparty/kfusion.h
	)
	
    cuda_add_library(${appname}-cuda    STATIC  ${kfusion_cuda_srcs} OPTIONS -O3 )
	target_link_libraries(${appname}-cuda   ${common_libraries})	
    add_version(${appname} cuda    "" "")

endif(CUDA_FOUND)

 # ----------------- PENCIL OpenCL VERSION -----------------

if (PENCIL_FOUND)
  message(STATUS "Adding a PENCIL version")
  pencil_compile(kfusion_kernel_srcs
    SOURCES src/pencil/pencil_kernels.c
    OPTIONS --no-private-memory --no-shared-memory --isl-schedule-max-coefficient=1 "--opencl-include-file=${CMAKE_CURRENT_SOURCE_DIR}/src/pencil/cl_kernel_vector.cl" --no-opencl-print-kernel-types
    )

  add_library(ppcgKernelsHost ${kfusion_kernel_srcs})
  target_include_directories(ppcgKernelsHost PRIVATE ${PENCIL_INCLUDE_DIRS})
  target_compile_definitions(ppcgKernelsHost PRIVATE ${PENCIL_DEFS})
  target_compile_options(ppcgKernelsHost PRIVATE -Wno-attributes)
  target_link_libraries(ppcgKernelsHost ${PENCIL_LIBS})
  set_target_properties(ppcgKernelsHost PROPERTIES C_STANDARD 99)
  set_target_properties(ppcgKernelsHost PROPERTIES C_STANDARD_REQUIRED ON)

  add_library(${appname}-pencilCL src/pencil/kernels.cpp)
  target_include_directories(${appname}-pencilCL PRIVATE ${PENCIL_INCLUDE_DIRS})
  target_link_libraries(${appname}-pencilCL ppcgKernelsHost)
  add_version(${appname} pencilCL "" "")
endif (PENCIL_FOUND)

#  ----------------- thirdparty ----------------- 

add_subdirectory( thirdparty    )


