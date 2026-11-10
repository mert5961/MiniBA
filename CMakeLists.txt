cmake_minimum_required(VERSION 3.20)
project(MiniBA CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Ceres REQUIRED)

add_executable(miniba src/miniba.cpp)
target_link_libraries(miniba PRIVATE Ceres::ceres)
target_compile_options(miniba PRIVATE -Wall -Wextra)

add_executable(test_synthetic tests/test_synthetic.cpp)
target_link_libraries(test_synthetic PRIVATE Ceres::ceres)

enable_testing()
add_test(NAME synthetic COMMAND test_synthetic)
