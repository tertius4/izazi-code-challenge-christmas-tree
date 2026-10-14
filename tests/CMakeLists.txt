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

add_executable(csv_tests
  args_test.cpp
  csv_reader_test.cpp
  csv_writer_test.cpp
)
target_link_libraries(csv_tests PRIVATE packer_lib GTest::gtest_main)
target_compile_options(csv_tests PRIVATE -Wall -Wextra)
target_compile_definitions(csv_tests PRIVATE
  TREE_CSV_PATH="${CMAKE_SOURCE_DIR}/data/input/tree.csv"
)

include(GoogleTest)
gtest_discover_tests(csv_tests)
