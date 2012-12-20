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
  acg_add_executable (small_miqp WIN32 ${sources} ${headers} )
else ()
  acg_add_executable (small_miqp ${sources} ${headers} )
endif ()

# enable rpath linking
set_target_properties(small_miqp PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

target_link_libraries (small_miqp
  CoMISo
  ${COMISO_LINK_LIBRARIES}
)