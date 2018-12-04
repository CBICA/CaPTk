function(copyToMacOSX target application)
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy \"${CMAKE_BINARY_DIR}/${application}\" 
        \"${CMAKE_BINARY_DIR}/${target}.app/Contents/Resources/bin/${application}\"
        COMMENT "Copying ${application}..."
    )
endfunction()