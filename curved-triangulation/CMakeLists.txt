
### OPENMESH
include_directories("D:/library/OpenMesh/OpenMesh9/src/")
link_directories("D:/library/OpenMesh/OpenMesh9/lib/")
add_definitions(-D_USE_MATH_DEFINES)

find_package(OpenMP REQUIRED)

aux_source_directory(. SOURCE_FILES)
add_executable(curved-triangulation ${SOURCE_FILES})

target_link_libraries(curved-triangulation edge-drawing)
target_link_libraries(curved-triangulation triangle)
target_link_libraries(curved-triangulation debug OpenMeshCored.lib optimized OpenMeshCore.lib)
target_link_libraries(curved-triangulation OpenMP::OpenMP_CXX)
