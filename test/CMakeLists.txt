set(test_names
    test_tv1_2d
    #other_test
    )

add_executable(test_tv1_2d "test_tv1_2d.cpp")
add_test(NAME tv1_2d COMMAND test_tv1_2d)

# To add more tests: append a name to test_names, and add_executable/add_test with that name
# Tests can be run with `ctest` in the build directory,
# it accepts regex, and different verbosity: ctest -R tv1 -VV
# add_executable(other_test "a_test_source_file.cpp")
# add_test(NAME a_name_that_makes_sense COMMAND other_test)

# Link proxTV and include_directory for each executable.
foreach(test_executable ${test_names})
    target_include_directories(${test_executable} PUBLIC
        $<BUILD_INTERFACE:${proxTV_BINARY_DIR}/src/include>)
    target_compile_definitions(${test_executable} PUBLIC NOMATLAB)
    target_link_libraries(${test_executable} PUBLIC proxTV)
endforeach()
