cmake_minimum_required(VERSION 2.8)

project(sph)

find_package(OpenGL REQUIRED)
find_package(OpenGL     REQUIRED)
find_package(GLUT       REQUIRED)

# Added fro openMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 11)

add_executable(sph ${CMAKE_CURRENT_SOURCE_DIR}/src/Visualizion_SPH.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/field2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/field3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/particlesystem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/wall.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/glvu.cpp
    )

# List of cpu sources
file(
    GLOB_RECURSE
    sources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)


file(GLOB_RECURSE app_sources *.cpp)

file(
    GLOB_RECURSE
    headers
    ${CMAKE_CURRENT_SOURCE_DIR}/include/*
)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_BUILD_TYPE}/${lib_dir})

target_link_libraries(sph ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )

