add_library( minimization OBJECT GslMinimizerHolder.cpp GslMinimizerHolder.hpp GslMinimizerRetriever.cpp GslMinimizerRetriever.hpp Minimizable.cpp Minimizable.hpp Minimizer.cpp Minimizer.hpp package.hpp )

add_library( TestMinimizable OBJECT TestMinimizable.cpp TestMinimizable.hpp )
foreach ( unittest GslMinimizerHolderTest MinimizableTest MinimizerTest )
	add_executable( ${unittest} ${unittest}.cpp $<TARGET_OBJECTS:TestSuite> $<TARGET_OBJECTS:TestMinimizable> $<TARGET_OBJECTS:minimization> $<TARGET_OBJECTS:linalg> )
	target_link_libraries( ${unittest} unit++ blas gsl )
	add_test( ${unittest} ${unittest} )
endforeach ( unittest )
