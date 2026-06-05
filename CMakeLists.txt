cmake_minimum_required(VERSION 3.14)
project(BarycentricOpenCL)

option(DISABLE_DEMO "Disable the demo program (which depends on tomlplusplus)" OFF)

include(FetchContent)

FetchContent_Declare(
    OpenCLWrapper
    GIT_REPOSITORY https://github.com/ProjectPhysX/OpenCL-Wrapper.git
    GIT_TAG        master
)

FetchContent_GetProperties(OpenCLWrapper)
if(NOT openclwrapper_POPULATED)
    FetchContent_MakeAvailable(OpenCLWrapper)
endif()

add_library(OpenCL_Wrapper INTERFACE)

target_include_directories(OpenCL_Wrapper SYSTEM INTERFACE ${openclwrapper_SOURCE_DIR}/src)

SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )

IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "-Wall -O2")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.45.0 CONFIG COMPONENTS graph) #graph filesystem iostreams program_options REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(fmt)
find_package(OpenCL REQUIRED)

include_directories(${CMAKE_CURRENT_BINARY_DIR} "${PROJECT_SOURCE_DIR}/src"
    ${Boost_INCLUDE_DIRS})

file(GLOB HEADER_FILES src/*.h)
file(GLOB CXX_FILES src/*.cpp)

add_library(${PROJECT_NAME} ${OS_BUNDLE} ${HEADER_FILES} ${CXX_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC OpenCL_Wrapper Eigen3::Eigen OpenCL::OpenCL ${Boost_LIBRARIES} fmt::fmt)

if(NOT DISABLE_DEMO)
    add_subdirectory(demo)
endif()
