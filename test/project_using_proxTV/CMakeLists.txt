cmake_minimum_required(VERSION 2.8.12)
project(dummy)
find_package(proxTV)
add_executable(dummy dummy.cpp)
# All the requirements of proxTV: compiler definitions, include_directories, linker flags (INTERFACE) are propagated to the target.
target_link_libraries(dummy PUBLIC proxTV::proxTV )
