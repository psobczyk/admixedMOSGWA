add_library( main OBJECT
	Exception.cpp
	Exception.hpp
	Firthimizer.cpp
	Firthimizer.hpp
	Individual.cpp
	Individual.hpp
	QRuncher.cpp
	QRuncher.hpp
	SNP.cpp
	SNP.hpp
)

add_subdirectory( io )
add_subdirectory( io/hdf5 )
add_subdirectory( linalg )
add_subdirectory( search )
add_subdirectory( search/egreedy )
add_subdirectory( search/memetica )
add_subdirectory( parser )
add_subdirectory( valuation )
add_subdirectory( util )
add_subdirectory( logging )
add_subdirectory( minimization )
add_subdirectory( lookup )

set( core_libraries
	main.cpp
	main_admixture.cpp
	Pool.cpp
	Pool.hpp
	buildinfo.cpp
	buildinfo.hpp
	GenotypeFreq.cpp
	GenotypeFreq.hpp
	Helpfull.cpp
	Helpfull.hpp
	MData.cpp
	MData.hpp
	Model.cpp
	Model.hpp
	Parameter.cpp
	Parameter.hpp
	PermSort.hpp
	SortVec.cpp
	SortVec.hpp
	ScoreTestShortcut.cpp
	ScoreTestShortcut.hpp
	AMdata.cpp
	AMdata.hpp
	AModel.cpp
	AModel.hpp
	$<TARGET_OBJECTS:main>
	$<TARGET_OBJECTS:io>
	$<TARGET_OBJECTS:linalg>
	$<TARGET_OBJECTS:search>
	$<TARGET_OBJECTS:egreedy>
	$<TARGET_OBJECTS:memetica>
	$<TARGET_OBJECTS:parser>
	$<TARGET_OBJECTS:valuation>
	$<TARGET_OBJECTS:util>
	$<TARGET_OBJECTS:logging>
	$<TARGET_OBJECTS:minimization>
	$<TARGET_OBJECTS:lookup>
)

find_package( BLAS REQUIRED )
find_package( HDF5 REQUIRED )
add_executable( MOSGWA ${core_libraries} $<TARGET_OBJECTS:io_hdf5> )
target_link_libraries( MOSGWA gsl gomp z pthread ${BLAS_LIBRARIES} ${HDF5_LIBRARIES} )
