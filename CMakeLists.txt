cmake_minimum_required(VERSION 3.5)

project(Matrix4D ${PROJECT_NAME} LANGUAGES CXX )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_INCLUDE_DIR ${_VCPKG_INSTALLED_DIR}/x64-windows-static/include)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost COMPONENTS  REQUIRED)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bigobj -ltbb -pthread")

 find_package(TBB CONFIG REQUIRED)

file(GLOB SOURCES
    include/*.hpp
    sources/*.cpp)
add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries (${PROJECT_NAME} ${Boost_LIBRARIES}  TBB::tbb)
