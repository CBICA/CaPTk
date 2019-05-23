find_package(Qt5Core REQUIRED)

# Retrieve the absolute path to qmake and then use that path to find
# the macdeployqt binary
get_target_property(_qmake_executable Qt5::qmake IMPORTED_LOCATION)
get_filename_component(_qt_bin_dir "${_qmake_executable}" DIRECTORY)
find_program(MACDEPLOYQT_EXECUTABLE macdeployqt HINTS "${_qt_bin_dir}")

# Add commands that copy the required Qt files to the application bundle
function(deployitksnap)
   execute_process(
        COMMAND sudo "${MACDEPLOYQT_EXECUTABLE}"
            "${CMAKE_BINARY_DIR}/../src/applications/individualApps/itksnap/ITK-SNAP.app" -always-overwrite
        COMMENT "Deploying Qt for ITK-SNAP..."
    )
endfunction()

mark_as_advanced(MACDEPLOYQT_EXECUTABLE)