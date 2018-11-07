set(LIB_INSTALL_DIR lib CACHE STRING "Install folder for lib (relative")
set(INCLUDE_INSTALL_DIR include CACHE STRING "Install folder for headers (relative")
set(EXECUTABLE_INSTALL_DIR bin CACHE STRING "Install folder for executables (relative")

set(headers
    condat_fast_tv.h
    general.h
    johnsonRyanTV.h
    lapackFunctionsWrap.h
    LPopt.h
    TVmacros.h
    TVopt.h
    utils.h
    )

set(sources
    condat_fast_tv.cpp
    johnsonRyanTV.cpp
    lapackFunctionsWrap.cpp
    LPopt.cpp
    TV2Dopt.cpp
    TV2DWopt.cpp
    TVgenopt.cpp
    TVL1opt.cpp
    TVL1opt_hybridtautstring.cpp
    TVL1opt_kolmogorov.cpp
    TVL1opt_tautstring.cpp
    TVL1Wopt.cpp
    TVL2opt.cpp
    TVLPopt.cpp
    TVNDopt.cpp
    utils.cpp
    )

add_library(proxTV ${sources})

target_compile_definitions(proxTV PUBLIC NOMATLAB)
if(proxTV_USE_LAPACK)
    target_compile_definitions(proxTV PUBLIC PROXTV_USE_LAPACK)
endif()

# Create directory with headers at build and install time.
file(COPY ${headers}
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/include")
install(FILES ${headers}
    DESTINATION ${INCLUDE_INSTALL_DIR})

target_include_directories(proxTV PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
    )

# Pthreads
target_link_libraries(proxTV PUBLIC
    ${CMAKE_THREAD_LIBS_INIT}
    )

# Lapacke, Lapack
if(proxTV_USE_LAPACK)
    target_include_directories(proxTV PUBLIC
    ${LAPACKE_INCLUDE_DIRS}
    )
    target_link_libraries(proxTV PUBLIC
        ${LAPACK_LIBRARIES}
        ${LAPACKE_LIBRARIES}
        )
else()
    target_link_libraries(proxTV PUBLIC Eigen3::Eigen)
endif()

# OpenMP
target_compile_options(proxTV PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(proxTV PUBLIC ${OpenMP_CXX_LIBRARIES})

install(TARGETS proxTV EXPORT proxTVTargets
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${EXECUTABLE_INSTALL_DIR}
    INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
    )

# This saves as to develop a hacky FindproxTV for other to use the library.
# A regular find_package(proxTV) will look for these targets.
# It needs an extra proxTVConfig.cmake to handle dependencies (provided in cmake directory)
install(EXPORT proxTVTargets
    FILE proxTVTargets.cmake
    NAMESPACE proxTV::
    DESTINATION ${LIB_INSTALL_DIR}/cmake/proxTV
    )

# Generate and install proxTVConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file("proxTVConfigVersion.cmake"
    VERSION ${proxTV_VERSION}
    COMPATIBILITY SameMajorVersion
    )
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/proxTVConfigVersion.cmake
    DESTINATION ${LIB_INSTALL_DIR}/cmake/proxTV
    )
