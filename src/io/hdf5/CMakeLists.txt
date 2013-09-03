add_library( io_hdf5 OBJECT Dataset.cpp Dataset.hpp Dataspace.hpp Datatype.cpp Datatype.hpp DoubleList.cpp DoubleList.hpp DoubleTable.cpp DoubleTable.hpp File.cpp File.hpp Id.cpp Id.hpp Object.hpp StringList.cpp StringList.hpp )
add_dependencies( io_hdf5 Dataspace.tpl Object.tpl )

add_library( Hdf5TestSuite OBJECT Hdf5TestSuite.cpp Hdf5TestSuite.hpp )
foreach ( unittest DoubleListTest DoubleTableTest ObjectTest StringListTest )
	add_executable( ${unittest} ${unittest}.cpp $<TARGET_OBJECTS:TestSuite> $<TARGET_OBJECTS:Hdf5TestSuite> $<TARGET_OBJECTS:main> $<TARGET_OBJECTS:io> $<TARGET_OBJECTS:io_hdf5> $<TARGET_OBJECTS:linalg> $<TARGET_OBJECTS:search> $<TARGET_OBJECTS:imputation> $<TARGET_OBJECTS:parser> $<TARGET_OBJECTS:valuation> $<TARGET_OBJECTS:util> $<TARGET_OBJECTS:minimization> $<TARGET_OBJECTS:lookup> )
	target_link_libraries( ${unittest} unit++ blas gsl hdf5 gomp z pthread )
	add_test( ${unittest} ${unittest} )
endforeach ( unittest )
