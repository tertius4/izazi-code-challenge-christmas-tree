cmake_minimum_required(VERSION 3.16)
project(packer CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(PACKER_BUILD_TESTS "Build unit tests" ON)

add_library(packer_lib
  src/io/csv_reader.cpp
  src/io/csv_writer.cpp
)
target_include_directories(packer_lib PUBLIC src)
target_compile_options(packer_lib PRIVATE -Wall -Wextra)

add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE packer_lib)
target_compile_options(main PRIVATE -Wall -Wextra)

if(PACKER_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
