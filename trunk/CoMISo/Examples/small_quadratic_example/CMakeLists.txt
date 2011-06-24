include (ACGCommon)
find_package(CoMISo)
find_package(SUITESPARSE)
find_package(MPI)

include_directories (
  ..
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${COMISO_INCLUDE_DIR}
  ${SUITESPARSE_INCLUDE_DIRS}
)

link_directories (
  ${SUITESPARSE_LIBRARY_DIRS}
  ${TAO_LIBRARY_DIR}
  ${PETSC_LIBRARY_DIR}
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
  acg_add_executable (small_quadratic_solver WIN32 ${sources} ${headers} )
elseif (APPLE)
  # generate bundle on mac
  acg_add_executable (small_quadratic_solver MACOSX_BUNDLE ${sources} ${headers} )
else ()
  acg_add_executable (small_quadratic_solver ${sources} ${headers} )
endif ()

# enable rpath linking
set_target_properties(small_quadratic_solver PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

target_link_libraries (small_quadratic_solver
  CoMISo
  ${SUITESPARSE_LIBRARIES}
  ${MPI_LIBRARY}
  ${PETSC_LIBRARY}
  ${TAO_LIBRARY}
)

if (APPLE)
   # create bundle in "Build" directory and set icon
   # no install needed here, because the whole bundle will be installed in the
   # toplevel CMakeLists.txt
   set_target_properties (
      Example PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build"
      MACOSX_BUNDLE_INFO_STRING "CoMISo small_quadratic_solver"
   )
endif ()


