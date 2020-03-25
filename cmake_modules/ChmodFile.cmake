function(chmodFile target application)
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND sudo chmod +x
    \"${CMAKE_BINARY_DIR}/${target}.app/Contents/Resources/bin/${application}\"
    COMMENT "Chmod +x ${application}"
  )
endfunction()