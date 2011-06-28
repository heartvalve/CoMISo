cmake_minimum_required (VERSION 2.6)

project(CoMISo)

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

include (ACGCommon)

acg_qt4 ()
# change to 0 if QT should not be used
set( WANT_COMISO_QT 1 )
if( QT4_FOUND)
  #message( WARNING " QT4 FOUND" )
  if( WANT_COMISO_QT )
    add_definitions (-DQT4_FOUND)
  #  message( WARNING " USING QT4" )
  endif ()
endif ()

acg_get_version ()

include (ACGOutput)

find_package (GMM)
if (GMM_FOUND)
  set (COMISO_GMM_CONFIG_FILE_SETTINGS "#define COMISO_GMM_AVAILABLE 1" )
else()
  set (COMISO_GMM_CONFIG_FILE_SETTINGS "#define COMISO_GMM_AVAILABLE 0" )
  message (FATAL_ERROR "GMM not found!")
endif ()

find_package (BLAS)
if (BLAS_FOUND )
  set (COMISO_BLAS_CONFIG_FILE_SETTINGS "#define COMISO_BLAS_AVAILABLE 1" )
else()
  set (COMISO_BLAS_CONFIG_FILE_SETTINGS "#define COMISO_BLAS_AVAILABLE 0" )
  message (FATAL_ERROR "BLAS not found!")
endif ()

#added by CAT
find_package (SUITESPARSE)
if (SUITESPARSE_FOUND )
  set (COMISO_SUITESPARSE_CONFIG_FILE_SETTINGS "#define COMISO_SUITESPARSE_AVAILABLE 1" )
else ()
  message (STATUS "SUITESPARSE not found!")
  set (COMISO_SUITESPARSE_CONFIG_FILE_SETTINGS "#define COMISO_SUITESPARSE_AVAILABLE 0" )
endif ()

find_package (MPI)
if (MPI_FOUND )
  set (COMISO_MPI_CONFIG_FILE_SETTINGS "#define COMISO_MPI_AVAILABLE 1" )
else ()
  message (STATUS "MPI not found!")
  set (COMISO_MPI_CONFIG_FILE_SETTINGS "#define COMISO_MPI_AVAILABLE 0" )
endif ()

find_package (PETSC)
if (PETSC_FOUND AND MPI_FOUND)
  set (COMISO_PETSC_CONFIG_FILE_SETTINGS "#define COMISO_PETSC_AVAILABLE 1" )
else ()
  message (STATUS "PETSC not found!")
  set (COMISO_PETSC_CONFIG_FILE_SETTINGS "#define COMISO_PETSC_AVAILABLE 0" )
endif ()


find_package (TAO)
if (TAO_FOUND AND PETSC_FOUND AND MPI_FOUND)
  set (COMISO_TAO_CONFIG_FILE_SETTINGS "#define COMISO_TAO_AVAILABLE 1" )
else ()
  message (STATUS "TAO not found!")
  set (COMISO_TAO_CONFIG_FILE_SETTINGS "#define COMISO_TAO_AVAILABLE 0" )
endif ()

find_package (IPOPT)
if (IPOPT_FOUND )
  set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 1" )
else ()
  message (STATUS "IPOPT not found!")
  set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 0" )
endif ()

if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/factored_solver/CMakeLists.txt" )
	add_subdirectory (Examples/factored_solver)
endif()

if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/quadratic_solver/CMakeLists.txt" )
	add_subdirectory (Examples/quadratic_solver)
endif()
if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/test2/CMakeLists.txt" )
	add_subdirectory (Examples/test2)
endif()
if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_quadratic_example/CMakeLists.txt" )
	add_subdirectory (Examples/small_quadratic_example)
endif()
if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_factored_example/CMakeLists.txt" )
	add_subdirectory (Examples/small_factored_example)
endif()

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${GMM_INCLUDE_DIR}
  ${SUITESPARSE_INCLUDE_DIRS}
  ${MPI_INCLUDE_PATH}
  ${PETSC_INCLUDE_DIRS}
  ${TAO_INCLUDE_DIRS}
  ${IPOPT_INCLUDE_DIR}
)

# generate dllexport macros on windows
if (WIN32)
  add_definitions(-DCOMISODLL)
endif ()


link_directories (
  ${SUITESPARSE_LIBRARY_DIRS}
#  ${MPI_LIBRARY_DIR}
  ${PETSC_LIBRARY_DIR}
  ${TAO_LIBRARY_DIR}
  ${IPOPT_LIBRARY_DIR}
)

# source code directories
set (directories 
  .
  Solver
  NSolver
  Config
  Utils
  QtWidgets
)

# generate dllexport macros on windows
if (WIN32)
  add_definitions(-DCOMISODLL)
endif ()

if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
  )
endif ()

# collect all header,source and ui files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})
acg_append_files (ui "*.ui" ${directories})


macro (of_list_filter _list)
  if (WIN32)
   foreach (_element ${${_list}})
     if (_element MATCHES "gnuplot_i\\.(cc|hh)$")
       list (REMOVE_ITEM ${_list} ${_element})
     endif ()
   endforeach ()
  endif ()
endmacro ()

of_list_filter  ( headers )
of_list_filter  ( sources )


# remove template cc files from source file list
acg_drop_templates (sources)

if( QT4_FOUND)
# genereate uic and moc targets
acg_qt4_autouic (uic_targets ${ui})
acg_qt4_automoc (moc_targets ${headers})
endif()
acg_add_library (CoMISo SHARED ${uic_targets} ${sources} ${headers} ${moc_targets})

if (NOT APPLE)
	target_link_libraries (CoMISo 
	  ${QT_LIBRARIES}
	  ${SUITESPARSE_LIBRARIES}
	  ${BLAS_LIBRARIES}
	)
else(NOT APPLE)
	target_link_libraries (CoMISo
	  ${QT_LIBRARIES}
	  ${SUITESPARSE_LIBRARIES}
	  ${BLAS_LIBRARIES}
	  ${PETSC_LIBRARIES}
	  ${TAO_LIBRARIES}
	  ${MPI_LIBRARIES}
	  ${IPOPT_LIBRARIES}
	)
endif(NOT APPLE)

# display results
acg_print_configure_header (COMISO "CoMISo")

# write config file
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/Config/config.hh.in"
     "${CMAKE_CURRENT_SOURCE_DIR}/Config/config.hh" @ONLY IMMEDIATE)

