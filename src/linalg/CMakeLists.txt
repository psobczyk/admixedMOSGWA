add_library( linalg OBJECT AutoMatrix.cpp AutoMatrix.hpp AutoPermutation.cpp AutoPermutation.hpp AutoVector.cpp AutoVector.hpp Matrix.cpp Matrix.hpp package.cpp package.hpp Permutation.cpp Permutation.hpp Vector.cpp Vector.hpp )

foreach ( unittest packageTest AutoPermutationTest AutoVectorTest AutoMatrixTest )
	add_executable( ${unittest} ${unittest}.cpp $<TARGET_OBJECTS:TestSuite> $<TARGET_OBJECTS:linalg> )
	target_link_libraries( ${unittest} unit++ blas gsl )
	add_test( ${unittest} ${unittest} )
endforeach ( unittest )
