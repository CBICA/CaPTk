function(changeDylib target application libname libpath)
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND 
    sudo install_name_tool -change \"${libpath}/${libname}\" \"@executable_path/../../Frameworks/${libname}\" \"${CMAKE_BINARY_DIR}/${application}\"

    COMMENT "Change ${libname} dylib path ${application}"
  )
endfunction()