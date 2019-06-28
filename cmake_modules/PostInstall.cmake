# Add symlink on osx and linux
IF (NOT WIN32)
  MESSAGE("Adding a symlink")
  INSTALL(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/bin/${EXE_NAME} {CMAKE_INSTALL_PREFIX}/bin/captk)")
ENDIF()
