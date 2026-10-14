cmake_minimum_required(VERSION 3.16)
project(packer CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE packer_lib)

add_library(packer_lib
  src/io/csv_reader.cpp
  src/io/csv_writer.cpp
)
target_include_directories(packer_lib PUBLIC src)
target_compile_options(packer_lib PRIVATE -Wall -Wextra)