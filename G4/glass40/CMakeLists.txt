cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(scint)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)

add_executable(scint ${sources} ${headers})
target_link_libraries(scint ${Geant4_LIBRARIES})

file(GLOB SCRIPTS ${PROJECT_SOURCE_DIR}/*.sh)
file(COPY ${SCRIPTS} DESTINATION ${PROJECT_BINARY_DIR})
file(GLOB PYSCRIPTS ${PROJECT_SOURCE_DIR}/*.py)
file(COPY ${PYSCRIPTS} DESTINATION ${PROJECT_BINARY_DIR})
