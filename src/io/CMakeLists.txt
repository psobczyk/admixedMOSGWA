add_library( io OBJECT Hdf5Constants.cpp Hdf5Constants.hpp Hdf5Input.cpp Hdf5Input.hpp Hdf5Output.cpp Hdf5Output.hpp InputAdapter.cpp InputAdapter.hpp Input.cpp Input.hpp Output.cpp Output.hpp PlinkInput.cpp PlinkInput.hpp )

foreach ( unittest PlinkInputTest Hdf5InputTest Hdf5OutputTest )
	add_executable( ${unittest} ${unittest}.cpp $<TARGET_OBJECTS:TestSuite> $<TARGET_OBJECTS:main> $<TARGET_OBJECTS:io> $<TARGET_OBJECTS:io_hdf5> $<TARGET_OBJECTS:linalg> $<TARGET_OBJECTS:search> $<TARGET_OBJECTS:imputation> $<TARGET_OBJECTS:parser> $<TARGET_OBJECTS:valuation> $<TARGET_OBJECTS:util> $<TARGET_OBJECTS:minimization> $<TARGET_OBJECTS:lookup> )
	target_link_libraries( ${unittest} unit++ blas gsl hdf5 gomp z pthread )
	add_test( ${unittest} ${unittest} )
endforeach ( unittest )
