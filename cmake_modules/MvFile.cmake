function(mvFile target app appname) 
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND sudo cp -R
    \"${CMAKE_BINARY_DIR}/externalApps/binaries_macos/${app}/${appname}.app\" \"${CMAKE_BINARY_DIR}/${target}.app/Contents/Resources/bin\"
    COMMENT "Copy ${app}"
  )
endfunction()
