cmake_minimum_required(VERSION 3.5)
project(emosqp C)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/public)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/private)

# Source files
file(GLOB OSQP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
file(GLOB WORKSPACE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*workspace.c)

# Static library
add_library(emosqpstatic STATIC ${OSQP_SRC} ${WORKSPACE_SRC})
target_link_libraries(emosqpstatic m)

# Example executable
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/emosqp.c)
    add_executable(emosqp ${CMAKE_CURRENT_SOURCE_DIR}/emosqp.c)
    target_link_libraries(emosqp emosqpstatic)
endif()
