add_executable (ibexsolve ibexsolve.cpp parse_args.h)
add_executable (ibexopt ibexopt.cpp parse_args.h)
add_executable (ibexminlpopt ibexminlpopt.cpp parse_args.h)
target_link_libraries (ibexsolve ibex)
target_link_libraries (ibexopt ibex)
target_link_libraries (ibexminlpopt ibex)
set_target_properties (ibexsolve PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                      "${CMAKE_BINARY_DIR}/bin")
set_target_properties (ibexopt PROPERTIES RUNTIME_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/bin")
set_target_properties (ibexminlpopt PROPERTIES RUNTIME_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/bin")

install (TARGETS ibexsolve DESTINATION ${CMAKE_INSTALL_BINDIR}
                           COMPONENT binaries)
install (TARGETS ibexopt DESTINATION ${CMAKE_INSTALL_BINDIR}
                           COMPONENT binaries)
install (TARGETS ibexminlpopt DESTINATION ${CMAKE_INSTALL_BINDIR}
                           COMPONENT binaries)
