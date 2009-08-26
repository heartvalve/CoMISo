include (ACGCommon)

include_directories (
  ..
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

link_directories (
/data/home1/zimmer/projects/ConstrainedSolver/package/CoMISo/build/Build/lib/CoMISo
)

# source code directories
set (directories 
  . 
)

# collect all header,source and ui files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})
acg_append_files (ui "*.ui" ${directories})

# remove template cc files from source file list
acg_drop_templates (sources)

# genereate uic and moc targets
acg_qt4_autouic (uic_targets ${ui})
acg_qt4_automoc (moc_targets ${headers})

if (WIN32)
  acg_add_executable (test2 WIN32 ${uic_targets} ${sources} ${headers} ${moc_targets})
elseif (APPLE)
  # generate bundle on mac
  acg_add_executable (test2 MACOSX_BUNDLE ${uic_targets} ${sources} ${headers} ${moc_targets})
else ()
  acg_add_executable (test2 ${uic_targets} ${sources} ${headers} ${moc_targets})
endif ()

target_link_libraries (test2
  CoMISo
  ${QT_LIBRARIES}
)

if (APPLE)
   # create bundle in "Build" directory and set icon
   # no install needed here, because the whole bundle will be installed in the
   # toplevel CMakeLists.txt
   set_target_properties (
      Exmaple PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build"
      MACOSX_BUNDLE_INFO_STRING "CoMISo test2"
   )
endif ()

