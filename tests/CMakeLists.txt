include_directories(${CMAKE_SOURCE_DIR}/others/)

find_package(Boost 1.46 COMPONENTS unit_test_framework random filesystem system timer)

add_executable(test_str test_str.cpp 
                        ${CMAKE_SOURCE_DIR}/others/string_utils.cpp)
target_link_libraries(test_str ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME Test_string_utils_module COMMAND test_str)


add_executable(test_expc test_expc.cpp)
target_link_libraries(test_expc ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} quadmath)
add_test(NAME Test_expc_function COMMAND test_expc)

add_executable(test_comb_points test_comb_points.cpp
                        ${CMAKE_SOURCE_DIR}/cryst_tools/comb_points.cpp)
target_link_libraries(test_comb_points ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_RANDOM_LIBRARY})
add_test(NAME Test_Comb_Points COMMAND test_comb_points)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


add_definitions(-DCT_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data_cryst_tools")



add_executable(test_cryst_tools test_cryst_tools.cpp
                        ${CMAKE_SOURCE_DIR}/cryst_tools/comb_points.cpp
                        ${CMAKE_SOURCE_DIR}/cryst_tools/cryst_tools.cpp)
target_link_libraries(test_cryst_tools ${Boost_LIBRARIES} PkgConfig::libOB)
add_test(NAME Test_Cryst_Tools COMMAND test_cryst_tools)


add_executable(test_combinatorics test_combinatorics.cpp)
target_link_libraries(test_combinatorics ${Boost_LIBRARIES})
add_test(NAME Test_Combinatorics COMMAND test_combinatorics)

