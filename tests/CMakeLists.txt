find_package(GTest QUIET)
if(NOT GTest_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
  )
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

add_executable(packer_tests
  args_test.cpp
  csv_reader_test.cpp
  csv_writer_test.cpp
  geometry_test.cpp
  io_test.cpp
  packing_test.cpp
  parallel_test.cpp
  tree_test.cpp
)
target_link_libraries(packer_tests PRIVATE packer_lib GTest::gtest_main)
target_compile_options(packer_tests PRIVATE -Wall -Wextra)
target_compile_definitions(packer_tests PRIVATE
  TREE_CSV_PATH="${CMAKE_SOURCE_DIR}/data/input/tree.csv"
)

include(GoogleTest)
gtest_discover_tests(packer_tests)
