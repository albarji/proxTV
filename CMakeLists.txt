cmake_minimum_required(VERSION 3.12)

project(proxTV
  VERSION 3.2.1
  DESCRIPTION "Toolbox for fast Total Variation proximity operators"
  )
message(STATUS "Configuring ${PROJECT_NAME}")
message(STATUS "  version: ${proxTV_VERSION}")
message(STATUS "  description: ${proxTV_DESCRIPTION}")

# Update CMake module path to lookup proxTV custom CMake modules
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#------------------------------------------------------------------------------
option(proxTV_INSTALL_DEVELOPMENT "Install development files" ON)
option(ENABLE_TESTING "Compile tests" OFF)
option(proxTV_USE_LAPACK "Use LAPACK and LAPACKE instead of EIGEN" ON)

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#------------------------------------------------------------------------------
# Install directories
if(NOT DEFINED proxTV_INSTALL_BIN_DIR)
   set(proxTV_INSTALL_BIN_DIR bin)
endif()
if(NOT DEFINED proxTV_INSTALL_LIB_DIR)
   set(proxTV_INSTALL_LIB_DIR lib)
endif()
if(NOT DEFINED proxTV_INSTALL_INCLUDE_DIR)
   set(proxTV_INSTALL_INCLUDE_DIR include)
endif()

#------------------------------------------------------------------------------
# External dependencies

if(proxTV_USE_LAPACK)
  find_package(LAPACKE REQUIRED)
  message(STATUS "Lapacke libraries: ${LAPACKE_LIBRARIES}")

  find_package(LAPACK REQUIRED)
  message(STATUS "Lapack libraries: ${LAPACK_LIBRARIES}")
else()
  find_package(Eigen3 REQUIRED)
  get_target_property(EIGEN_INCLUDE_DIR Eigen3::Eigen
    INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Eigen Found: ${EIGEN_INCLUDE_DIR}")
endif()

set(THREADS_PREFER_PTHREAD_FLAG 1)
find_package(Threads)

find_package(OpenMP)
message(STATUS "OpenMP found: ${OpenMP_FOUND}")

#------------------------------------------------------------------------------
# Add libraries

add_subdirectory(src)

#------------------------------------------------------------------------------
# Testing
if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

#------------------------------------------------------------------------------
# Configure proxTVConfigVersion.cmake common to build and install tree
include(CMakePackageConfigHelpers)
set(config_version_file "${proxTV_BINARY_DIR}/proxTVConfigVersion.cmake")
write_basic_package_version_file(${config_version_file}
    VERSION ${proxTV_VERSION}
    COMPATIBILITY SameMajorVersion
    )

#------------------------------------------------------------------------------
# Export 'proxTVTargets.cmake' for a build tree
export(TARGETS proxTV
  FILE ${PROJECT_BINARY_DIR}/proxTVTargets.cmake
  NAMESPACE proxTV::
  )

# Configure 'proxTVConfig.cmake' for a build tree
include(CMakePackageConfigHelpers)
set(build_config ${PROJECT_BINARY_DIR}/proxTVConfig.cmake)
configure_package_config_file(
    cmake/proxTVConfig.cmake.in
    ${build_config}
    INSTALL_DESTINATION ${PROJECT_BINARY_DIR}
    INSTALL_PREFIX ${PROJECT_BINARY_DIR}
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

#------------------------------------------------------------------------------
# Configure 'proxTVConfig.cmake' for an install tree
if(proxTV_INSTALL_DEVELOPMENT)
    set(install_config ${PROJECT_BINARY_DIR}/CMakeFiles/proxTVConfig.cmake)
    configure_package_config_file(
        cmake/proxTVConfig.cmake.in
        ${install_config}
        INSTALL_DESTINATION ${proxTV_INSTALL_LIB_DIR}/cmake/proxTV
        NO_CHECK_REQUIRED_COMPONENTS_MACRO
        )

    # Install 'proxTVTargets.cmake'
    install(EXPORT proxTVTargets
        FILE proxTVTargets.cmake
        NAMESPACE proxTV::
        DESTINATION ${proxTV_INSTALL_LIB_DIR}/cmake/proxTV
        COMPONENT Development
        )

    # Install config files
    install(
        FILES ${config_version_file} ${install_config}
        DESTINATION ${proxTV_INSTALL_LIB_DIR}/cmake/proxTV
        COMPONENT Development
        )
endif()
