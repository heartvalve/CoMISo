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
if (NOT GMM_FOUND)
  message (FATAL_ERROR "GMM not found!")
endif ()

find_package (Boost)
if (NOT Boost_FOUND)
  message (FATAL_ERROR "Boost not found!")
endif ()

find_package (BLAS)
if (NOT BLAS_FOUND )
  message (FATAL_ERROR "BLAS not found!")
endif ()

#added by CAT
find_package (SUITESPARSE)
if (NOT SUITESPARSE_FOUND )
  message (FATAL_ERROR "SUITESPARSE not found!")
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
  ${Boost_INCLUDE_DIR}
  ${GMM_INCLUDE_DIR}
  ${SUITESPARSE_INCLUDE_DIRS}
)

# generate dllexport macros on windows
if (WIN32)
  add_definitions(-DCOMISODLL)
endif ()


link_directories (
  ${SUITESPARSE_LIBRARY_DIRS}
)

# source code directories
set (directories 
  .
  Solver 
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

message( WARNING ${ui})

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

# genereate uic and moc targets
acg_qt4_autouic (uic_targets ${ui})
acg_qt4_automoc (moc_targets ${headers})

message( WARNING ${uic_targets})
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
	)
endif(NOT APPLE)

# display results
acg_print_configure_header (COMISO "CoMISo")

