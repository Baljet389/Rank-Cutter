

add_executable(rac_test testMain.cpp 
				"testArgs.h"
				"testArgs.cpp"
				"testSVD.h"
				"testSVD.cpp")

target_include_directories(rac_test PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../inc
    ${CMAKE_CURRENT_SOURCE_DIR}/../lib
)

if(MSVC)
	target_compile_options(rac_test
	PRIVATE /W4 
			/permissive- 
			$<$<CONFIG:Release>:/O2 
			/arch:AVX2 
			/fp:fast>
			)
else()
	target_compile_options(rac_test
	PRIVATE -Wall         
			-Wextra        
			-Wpedantic
			-march=native 
			-ffast-math 
			-O3 
	)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(rac_test PUBLIC OpenMP::OpenMP_CXX)
	endif()
endif()