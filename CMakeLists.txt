cmake_minimum_required(VERSION 3.5)

project(spatial_group ${PROJECT_NAME} LANGUAGES CXX )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(Boost_USE_MULTITHREADED ON)



find_package(Boost ${Boost_VERSION_MAJOR} REQUIRED timer#[[ regex unit_test_framework]]  )


find_package(TBB REQUIRED tbb)


add_compile_definitions(__TBB_DYNAMIC_LOAD_ENABLED=0)

file(GLOB SOURCES
    include/*.hpp
    sources/main.cpp)

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ltbb -fopenmp -pedantic -pthread")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=2900")


add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries (${PROJECT_NAME} ${Boost_LIBRARIES} TBB::tbb)

#set (TESTS_SOURCES sources/unit_test.cpp)
#set (TEST test_${PROJECT})
#add_executable (${TEST} ${TESTS_SOURCES})
#target_link_libraries (${TEST} ${PROJECT} ${Boost_LIBRARIES} TBB::tbb MPI::MPI_CXX)
##enable_testing ()
#add_test (${TEST} ${TEST})

