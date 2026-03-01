cmake_minimum_required(VERSION 3.16)

project(rac VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP)
include_directories(inc)

add_executable(${PROJECT_NAME}
	src/main.cpp
	src/stb_impl.cpp

	"inc/matrix.h"
	"inc/svd.h"
	
	"inc/qr.h"
	"inc/image.h"
	"inc/utility.h" 
	"inc/bidiagonal.h" 
	"inc/householder.h" 
	"inc/blockUpdate.h"
	"inc/args.h" 
	"src/args.cpp")

target_include_directories(${PROJECT_NAME} SYSTEM BEFORE PRIVATE lib)
add_compile_definitions(PROJECT_NAME=${PROJECT_NAME})
add_compile_definitions(PROJECT_VERSION=${PROJECT_VERSION})

if(MSVC)
	target_compile_options(${PROJECT_NAME}
	PRIVATE /W4 
			/permissive- 
			$<$<CONFIG:Release>:/O2 
			/arch:AVX2 
			/fp:fast>
			)
else()
	target_compile_options(${PROJECT_NAME}
	PRIVATE -Wall         
			-Wextra        
			-Wpedantic
			-march=native 
			-ffast-math 
			-O3 
			#-fsanitize=address -g
			#-Rpass-analysis=loop-vectorize
	)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
	endif()
	# target_link_options(${PROJECT_NAME} PRIVATE -fsanitize=address)
endif()

add_subdirectory(test)
  
