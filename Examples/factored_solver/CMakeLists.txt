include (ACGCommon)

include_directories (
  ..
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

link_directories (
 # ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/CoMISo
#/data/home1/zimmer/projects/ConstrainedSolver/package/CoMISo/build/Build/lib/CoMISo
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
  acg_add_executable (factored_solver WIN32 ${sources} ${headers} )
elseif (APPLE)
  # generate bundle on mac
  acg_add_executable (factored_solver MACOSX_BUNDLE ${sources} ${headers} )
else ()
  acg_add_executable (factored_solver ${sources} ${headers} )
endif ()

target_link_libraries (factored_solver
  CoMISo
)

if (APPLE)
   # create bundle in "Build" directory and set icon
   # no install needed here, because the whole bundle will be installed in the
   # toplevel CMakeLists.txt
   set_target_properties (
      Exmaple PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build"
      MACOSX_BUNDLE_INFO_STRING "CoMISo factored_solver"
   )
endif ()

