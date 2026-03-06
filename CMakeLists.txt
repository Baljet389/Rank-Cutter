cmake_minimum_required(VERSION 3.16)

project(rac VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP)

add_library(MathSettings INTERFACE)
add_compile_definitions(PROJECT_NAME=${PROJECT_NAME})
add_compile_definitions(PROJECT_VERSION=${PROJECT_VERSION})

target_include_directories(MathSettings INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
)
target_include_directories(MathSettings SYSTEM INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/lib
)

if(MSVC)
    target_compile_options(MathSettings INTERFACE
        /W4 /permissive- 
        $<$<CONFIG:Release>:/O2 /arch:AVX2 /fp:fast>
    )
else()
    target_compile_options(MathSettings INTERFACE
        -Wall -Wextra -Wpedantic
        -march=native -ffast-math -O3
    )
    if(OpenMP_CXX_FOUND)
        target_link_libraries(MathSettings INTERFACE OpenMP::OpenMP_CXX)
    endif()
endif()

add_executable(${PROJECT_NAME} src/main.cpp src/stb_impl.cpp src/args.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE MathSettings)


add_subdirectory(test)
  
