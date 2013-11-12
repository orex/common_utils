include_directories(${CMAKE_SOURCE_DIR}/others/)

add_executable(test_str test_str.cpp 
                        ${CMAKE_SOURCE_DIR}/others/string_utils.cpp)
target_link_libraries(test_str ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME Test_string_utils_module COMMAND test_str)


add_executable(test_expc test_expc.cpp)
target_link_libraries(test_expc ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} quadmath)
add_test(NAME Test_expc_function COMMAND test_expc)

find_package(Boost 1.46 COMPONENTS random)

add_executable(test_rnd test_rnd.cpp
                        ${CMAKE_SOURCE_DIR}/others/rnd_utils.cpp)
target_link_libraries(test_rnd ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_RANDOM_LIBRARY})
add_test(NAME Test_random_array COMMAND test_rnd)

