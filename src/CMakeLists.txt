add_library(
  pt STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/vec.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/random_state.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ray.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/sphere.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hit_record.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/camera.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp)
target_include_directories(pt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pt PRIVATE project_options project_warnings)

add_executable(${CMAKE_PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(
  ${CMAKE_PROJECT_NAME}
  PRIVATE project_options
          project_warnings
          fmt::fmt
          cpp-taskflow::cpp-taskflow
          pt)
