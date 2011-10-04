include (ACGCommon)
find_package(CoMISo)
find_package(SUITESPARSE)
find_package(MPI)
find_package (MUMPS)
find_package(IPOPT)

if (MUMPS_FOUND )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${MUMPS_INCLUDE_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES ${MUMPS_LIBRARY} )
endif ()

find_package (IPOPT)
if (IPOPT_FOUND AND MUMPS_FOUND)
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${IPOPT_INCLUDE_DIR} )
  list( APPEND COMISO_LINK_DIRECTORIES ${IPOPT_LIBRARY_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES ${IPOPT_LIBRARY} )
endif ()

include_directories (
  ..
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${COMISO_INCLUDE_DIR}
  ${SUITESPARSE_INCLUDE_DIRS}
  ${COMISO_INCLUDE_DIRECTORIES}
)

link_directories (
  ${SUITESPARSE_LIBRARY_DIRS}
  ${TAO_LIBRARY_DIR}
  ${PETSC_LIBRARY_DIR}
  ${COMISO_LINK_DIRECTORIES}
)

# source code directories
set (directories 
  . 
)

# collect all header and source files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})

# remove template cc files from source file list
acg_drop_templates (sources)

if (WIN32)
  acg_add_executable (small_factored_solver WIN32 ${sources} ${headers} )
elseif (APPLE)
  # generate bundle on mac
  acg_add_executable (small_factored_solver MACOSX_BUNDLE ${sources} ${headers} )
else ()
  acg_add_executable (small_factored_solver ${sources} ${headers} )
endif ()

# enable rpath linking
set_target_properties(small_factored_solver PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

target_link_libraries (small_factored_solver
  CoMISo
  ${SUITESPARSE_LIBRARIES}
  ${MPI_LIBRARY}
  ${PETSC_LIBRARY}
  ${TAO_LIBRARY}
  ${COMISO_LINK_LIBRARIES}
)

if (APPLE)
   # create bundle in "Build" directory and set icon
   # no install needed here, because the whole bundle will be installed in the
   # toplevel CMakeLists.txt
   set_target_properties (
      Example PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build"
      MACOSX_BUNDLE_INFO_STRING "CoMISo small_factored_solver"
   )
endif ()
