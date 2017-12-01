cmake_minimum_required(VERSION 2.8.9)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project (project8_data_gathering)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
set(SRC_PATH "Sources")

#------------------------------------------------------------
# Use C++ 11
#------------------------------------------------------------

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#------------------------------------------------------------
# Source Files
#------------------------------------------------------------
add_subdirectory(src)

set(SRC 
    ${project8_data_gathering_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PixelIterator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PixelIterator.h
    ) 



#------------------------------------------------------------
# External Dependencies
#------------------------------------------------------------
FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIR})

# For IDE folder hierarchy generation
set(SRC_PATH ".")
file(
    GLOB_RECURSE SRC_LIST 
    LIST_DIRECTORIES false
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${SRC_PATH}/*.c*"
    "${SRC_PATH}/*.h*"
    "${SRC_PATH}/*.txt*"
    "${SRC_PATH}/*.frag*"
    "${SRC_PATH}/*.vert*"
)

# Definitions 
add_definitions(-DSCENE_PATH="${CMAKE_SOURCE_DIR}/../Scenes/")
add_definitions(-DMODEL_PATH="${CMAKE_SOURCE_DIR}/../Models/")
add_definitions(-DTEXTURE_PATH="${CMAKE_SOURCE_DIR}/../Textures/")

# Main Target
add_executable(project8_data_gathering "${SRC}")
target_link_libraries(project8_data_gathering ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})
set_target_properties(project8_data_gathering PROPERTIES LINKER_LANGUAGE CXX)

# For IDE folder hierarchy generation
foreach(_source IN ITEMS ${SRC_LIST})
	get_filename_component(_source_path "${_source}" PATH)
	string(REPLACE "/" "\\" _source_path_msvc "${_source_path}")
	source_group("${_source_path_msvc}" FILES "${_source}")
endforeach() 
