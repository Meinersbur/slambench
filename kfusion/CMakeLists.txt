cmake_minimum_required(VERSION 2.8)

set(appname kfusion)

find_package(TooN   REQUIRED)
find_package(OpenCL )
find_package(CUDA   )
find_package(GLUT   )
find_package(OpenGL )
find_package(OpenNI2)

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

add_library(${appname}-cpp  src/cpp/kernels.cpp)
target_link_libraries(${appname}-cpp   ${common_libraries})	
add_version(${appname} cpp "" "")

 # ----------------- OPENMP VERSION ----------------- 

add_library(${appname}-openmp  src/cpp/kernels.cpp)
target_link_libraries(${appname}-openmp   ${common_libraries})	
SET_TARGET_PROPERTIES(${appname}-openmp PROPERTIES COMPILE_FLAGS "-fopenmp")
add_version(${appname} openmp "-fopenmp" "-fopenmp")


 #  ----------------- OCL VERSION ----------------- 
 
if (OPENCL_FOUND) 
   include_directories(${OPENCL_INCLUDE_DIRS})
    add_library(${appname}-opencl  src/opencl/kernels.cpp src/opencl/common_opencl.cpp)
    target_link_libraries(${appname}-opencl   ${common_libraries} ${OPENCL_LIBRARIES})	
    add_version(${appname} opencl "" "")
endif(OPENCL_FOUND)


 #  ----------------- CUDA VERSION ----------------- 
if (CUDA_FOUND) 


set(kfusion_cuda_srcs
	src/cuda/kernels.cu
	thirdparty/kfusion.h
	)
	
    cuda_add_library(${appname}-cuda    STATIC  ${kfusion_cuda_srcs} OPTIONS -O3 )
	target_link_libraries(${appname}-cuda   ${common_libraries})	
    add_version(${appname} cuda    "" "")

endif(CUDA_FOUND)

 # ----------------- PENCIL OpenCL VERSION -----------------

set(KFUSION_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(KFUSION_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(PENCILCC_EXECUTABLE CACHE FILEPATH "Path to pencilcc")

if (OPENCL_FOUND AND PENCILCC_EXECUTABLE)
	execute_process(COMMAND "${PENCILCC_EXECUTABLE}" --noselfupdate --show-cc-args=CPPFLAGS
		OUTPUT_VARIABLE PENCIL_CC_CPPFLAGS
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND "${PENCILCC_EXECUTABLE}" --noselfupdate --show-ld-args=LIBRARY_DIRS
		OUTPUT_VARIABLE PENCIL_LD_LIBRARY_DIRS
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND "${PENCILCC_EXECUTABLE}" --noselfupdate --show-ld-args=LIBS
		OUTPUT_VARIABLE PENCIL_LD_LIBS
		OUTPUT_STRIP_TRAILING_WHITESPACE)
#	execute_process(COMMAND "llvm-config" --libs
	execute_process(COMMAND "${PENCILCC_EXECUTABLE}" --noselfupdate --show-cc-args=INCLUDE_DIRS
		OUTPUT_VARIABLE PENCIL_CC_INCLUDE_DIRS
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	string(REPLACE " " ";" PENCIL_CC_INCLUDE_DIRS ${PENCIL_CC_INCLUDE_DIRS})

	add_custom_command(
		OUTPUT "${KFUSION_BINARY_DIR}/pencil_kernels_opt_host.c"
		DEPENDS "${KFUSION_SOURCE_DIR}/src/pencil/pencil_kernels.c"
		COMMAND "${PENCILCC_EXECUTABLE}" "${KFUSION_SOURCE_DIR}/src/pencil/pencil_kernels.c" --pencil-mode=ppcg --no-private-memory --no-shared-memory --isl-schedule-max-coefficient=1 -o "${KFUSION_BINARY_DIR}/pencil_kernels_opt_host.c" --show-commands --opencl-include-file=${CMAKE_CURRENT_SOURCE_DIR}/src/pencil/cl_kernel_vector.cl --no-opencl-print-kernel-types )

	#set(prl_path $ENV{PRL_PATH})
	#set(prl_src_path $ENV{PRL_SRC_PATH})
	#set(util_path $ENV{PENCIL_UTIL_HOME})
	#include_directories(${PENCIL_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} "${util_path}/include" "${prl_src_path}/include")
	#add_definitions(${PENCIL_CC_CPPFLAGS})

	add_library(ppcgKernelsHost "${KFUSION_BINARY_DIR}/pencil_kernels_opt_host.c")
	SET_TARGET_PROPERTIES(ppcgKernelsHost PROPERTIES COMPILE_FLAGS "-std=c99 ${PENCIL_CC_CPPFLAGS}")

	#add_library(ppcgKernelsNative "src/pencil/kernel_vector.c")
	#SET_TARGET_PROPERTIES(ppcgKernelsNative PROPERTIES COMPILE_FLAGS "-std=c99")
	#foreach (_dir IN LISTS PENCIL_CC_INCLUDE_DIRS)
	#	target_include_directories(ppcgKernelsNative PRIVATE ${_dir})
	#endforeach ()

	find_library(PRL_LIBRARY NAMES "${PENCIL_LD_LIBS}" HINTS ${PENCIL_LD_LIBRARY_DIRS})

	add_library(${appname}-pencilCL src/pencil/kernels.cpp)
	SET_TARGET_PROPERTIES(${appname}-pencilCL PROPERTIES COMPILE_FLAGS "${PENCIL_CC_CPPFLAGS}")
	target_link_libraries(${appname}-pencilCL ppcgKernelsHost ${PRL_LIBRARY} ${OPENCL_LIBRARIES} ${common_libraries})
	add_version(${appname} pencilCL "${PENCIL_CC_CPPFLAGS}" "${PRL_LIBRARY}")
endif()

#  ----------------- thirdparty ----------------- 

add_subdirectory( thirdparty    )


