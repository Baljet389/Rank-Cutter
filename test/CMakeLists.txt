
add_executable(rac_test testMain.cpp testArgs.cpp testSVD.cpp)
target_link_libraries(rac_test PRIVATE MathSettings)