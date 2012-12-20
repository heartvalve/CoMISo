include (ACGCommon)
include (CoMISoExample)

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
  acg_add_executable (quadratic_solver WIN32 ${sources} ${headers} )
else ()
  acg_add_executable (quadratic_solver ${sources} ${headers} )
endif ()

# enable rpath linking
set_target_properties(quadratic_solver PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

target_link_libraries (quadratic_solver
  CoMISo
  ${COMISO_LINK_LIBRARIES}
)
