cmake_minimum_required (VERSION 3.2)

project(common_utils)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-sign-compare -Wextra " )

find_package(Boost 1.46 COMPONENTS unit_test_framework)

if(Boost_UNIT_TEST_FRAMEWORK_FOUND)
  enable_testing()
else()
  message(STATUS "No Boost::Test found. All tests have been disabled.")
endif()  

#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
find_package(PkgConfig)
pkg_check_modules(libOB REQUIRED IMPORTED_TARGET openbabel-3)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(Boost_UNIT_TEST_FRAMEWORK_FOUND)
  add_subdirectory(tests)
endif()    
