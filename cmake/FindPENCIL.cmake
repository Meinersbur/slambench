 
include(CMakeParseArguments)
include(FindPackageHandleStandardArgs)


find_program(PENCILCC_EXECUTABLE pencilcc)

if (PENCILCC_EXECUTABLE)
  execute_process(COMMAND "${PENCILCC_EXECUTABLE}" "--show-listsep=\;" --show-cc-args=INCLUDE_DIRS
    OUTPUT_VARIABLE PENCIL_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  set(PENCIL_INCLUDE_DIRS ${PENCIL_INCLUDE_DIRS})
  execute_process(COMMAND "${PENCILCC_EXECUTABLE}" "--show-listsep=\;" --show-cc-args=DEFFLAGS
    OUTPUT_VARIABLE PENCIL_DEFS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  set(PENCIL_DEFS ${PENCIL_DEFS})
  execute_process(COMMAND "${PENCILCC_EXECUTABLE}" "--show-listsep=\;" --show-cc-args=LIBS
    OUTPUT_VARIABLE PENCIL_REL_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  set(PENCIL_REL_LIBS ${PENCIL_REL_LIBS})
  execute_process(COMMAND "${PENCILCC_EXECUTABLE}" "--show-listsep=\;" --show-cc-args=LIBRARY_DIRS
    OUTPUT_VARIABLE PENCIL_LIBRARY_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  set(PENCIL_LIBRARY_DIRS ${PENCIL_LIBRARY_DIRS})

  set(PENCIL_LIBS)
  foreach (lib IN LISTS PENCIL_REL_LIBS)
    find_library("PENCIL_LIB_${lib}" NAMES "${lib}" HINTS ${PENCIL_LIBRARY_DIRS})
    list(APPEND PENCIL_LIBS "${PENCIL_LIB_${lib}}" )
  endforeach ()
endif ()

find_package_handle_standard_args(PENCIL
  REQUIRED_VARS PENCILCC_EXECUTABLE
  )

function (pencil_compile var)
  cmake_parse_arguments(ARG "" "" "OPTIONS;SOURCES" ${ARGN})

  set(results)
  foreach(infile IN LISTS ARG_SOURCES)
    get_filename_component(infile_abspath "${infile}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    get_filename_component(basename "${infile}" NAME_WE)
    set(outfile "${basename}_host.c")
    get_filename_component(outfile_abspath "${outfile}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_BINARY_DIR}")

    add_custom_command(
      OUTPUT "${outfile_abspath}"
      DEPENDS "${infile_abspath}"
      COMMAND "${PENCILCC_EXECUTABLE}" "${infile_abspath}" --no-cc -c ${ARG_OPTIONS} -o "${outfile_abspath}"
      COMMENT "ppcg ${infile}"
      )
    list(APPEND results ${outfile_abspath})
  endforeach ()
  set("${var}" ${results} PARENT_SCOPE)
endfunction ()



if (OFF)
set(GMP_DEFINITIONS)

find_path(GMP_INCLUDE_DIR gmp.h
  PATHS ENV GMP_INCLUDE ENV GMP_DIR
  NO_DEFAULT_PATH
)
find_path(GMP_INCLUDE_DIR gmp.h)
mark_as_advanced(GMP_INCLUDE_DIR)
set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})

find_library(GMP_LIBRARY NAMES gmp mpir
  HINTS ENV GMP_LIB ENV GMP_DIR
  NO_DEFAULT_PATH
)
find_library(GMP_LIBRARY NAMES gmp mpir)
mark_as_advanced(GMP_LIBRARY)
set(GMP_LIBRARIES ${GMP_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmp DEFAULT_MSG GMP_LIBRARY GMP_INCLUDE_DIR)



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

  add_custom_command(
    OUTPUT "${KFUSION_BINARY_DIR}/pencil_kernels_opt_host.c"
    DEPENDS "${KFUSION_SOURCE_DIR}/src/pencil/pencil_kernels.c"
    COMMAND "${PENCILCC_EXECUTABLE}" "${KFUSION_SOURCE_DIR}/src/pencil/pencil_kernels.c" --pencil-mode=ppcg --no-private-memory --no-shared-memory --isl-schedule-max-coefficient=1 -o "${KFUSION_BINARY_DIR}/pencil_kernels_opt_host.c" --show-commands --opencl-include-file=${CMAKE_CURRENT_SOURCE_DIR}/src/pencil/cl_kernel_vector.cl --no-opencl-print-kernel-types )


	add_library(${appname}-pencilCL src/pencil/kernels.cpp)
	SET_TARGET_PROPERTIES(${appname}-pencilCL PROPERTIES COMPILE_FLAGS "${PENCIL_CC_CPPFLAGS}")
	target_link_libraries(${appname}-pencilCL ppcgKernelsHost ${PRL_LIBRARY} ${OPENCL_LIBRARIES} ${common_libraries})
	add_version(${appname} pencilCL "${PENCIL_CC_CPPFLAGS}" "${PRL_LIBRARY}")
endif()
endif()