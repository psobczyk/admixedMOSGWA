add_library( lookup OBJECT ModelIndex.cpp ModelIndex.hpp MutableModelIndex.cpp MutableModelIndex.hpp package.cpp package.hpp ResultStore.cpp ResultStore.hpp )

foreach ( unittest ModelIndexTest )
	add_executable( ${unittest} ${unittest}.cpp $<TARGET_OBJECTS:TestSuite> $<TARGET_OBJECTS:lookup> )
	target_link_libraries( ${unittest} unit++ blas gsl )
	add_test( ${unittest} ${unittest} )
endforeach ( unittest )
