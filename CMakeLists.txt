cmake_minimum_required (VERSION 3.1)

project("image-curved-triangulation")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_definitions(-D_CRT_SECURE_NO_WARNINGS)

add_subdirectory(extern/triangle)
#add_subdirectory(edge-drawing)

#add_subdirectory(edge-detection)
add_subdirectory(curved-triangulation)



