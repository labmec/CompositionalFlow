cmake_minimum_required(VERSION 3.14)

project (DivFreeBubbles LANGUAGES C CXX)


set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)


# Making pz "seeable" in Xcode
file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})


# creating the projects
add_subdirectory(targets)
