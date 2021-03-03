function(enable_openmp project_name)
  find_package(OpenMP)

  if(OpenMP_FOUND)
    target_link_libraries(${project_name} INTERFACE ${OpenMP_CXX_FLAGS})
  else()
    message("Could not find OpenMP on the system!")
  endif()
endfunction()
