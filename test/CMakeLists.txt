set(test_names
    test_tv1_2d
    #other_test
    )

add_executable(test_tv1_2d "test_tv1_2d.cpp")
add_test(NAME tv1_2d COMMAND test_tv1_2d)

if(proxTV_USE_LAPACK)
    find_package(Eigen3)
    if (TARGET Eigen3::Eigen)
        list(APPEND test_names test_use_eigen)
        add_executable(test_use_eigen "test_use_eigen.cpp")
        add_test(NAME use_eigen COMMAND test_use_eigen)
        target_link_libraries(test_use_eigen PUBLIC Eigen3::Eigen)
    else()
        message(WARNING "Eigen3 not found, test comparing lapack and eigen is disabled.")
    endif()
endif()
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
