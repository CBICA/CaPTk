function(mvFile target old new) 
    add_custom_command(TARGET ${target} POST_BUILD
    COMMAND mv
    \"${CMAKE_BINARY_DIR}/${target}.app/Contents/Resources/bin/${old}\" \"${CMAKE_BINARY_DIR}/${target}.app/Contents/Resources/bin/${new}\"
    COMMENT "Rename ${old} to ${new}"
    )
endfunction()