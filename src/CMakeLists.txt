add_library( main OBJECT Exception.cpp Exception.hpp Firthimizer.cpp Firthimizer.hpp Individual.cpp Individual.hpp QRuncher.cpp QRuncher.hpp SNP.cpp SNP.hpp types.hpp )

add_subdirectory( io )
add_subdirectory( io/hdf5 )
add_subdirectory( linalg )
add_subdirectory( search )
add_subdirectory( imputation )
add_subdirectory( parser )
add_subdirectory( valuation )
add_subdirectory( util )
add_subdirectory( minimization )
add_subdirectory( lookup )

add_custom_command(
  OUTPUT timestamp.txt
  COMMAND cmake -E echo_append "#define TIMESTAMP \"" > timestamp.txt
  COMMAND date "+%F %T %Z\"" >> timestamp.txt
  VERBATIM
)
add_library( buildinfo OBJECT buildinfo.cpp buildinfo.hpp timestamp.txt )

set( core_libraries main.cpp GenotypeFreq.cpp GenotypeFreq.hpp Helpfull.cpp Helpfull.hpp Log.cpp Log.hpp MData.cpp MData.hpp Model.cpp Model.hpp Parameter.cpp Parameter.hpp PermSort.hpp SortVec.cpp SortVec.hpp ScoreTestShortcut.cpp ScoreTestShortcut.hpp $<TARGET_OBJECTS:main> $<TARGET_OBJECTS:buildinfo> $<TARGET_OBJECTS:io> $<TARGET_OBJECTS:linalg> $<TARGET_OBJECTS:search> $<TARGET_OBJECTS:imputation> $<TARGET_OBJECTS:parser> $<TARGET_OBJECTS:valuation> $<TARGET_OBJECTS:util> $<TARGET_OBJECTS:minimization> $<TARGET_OBJECTS:lookup> )

set ( external_core_libraries blas gsl gomp z pthread )

find_package( BLAS REQUIRED )
find_package( HDF5 REQUIRED )
if ( HDF5_FOUND )
	message( STATUS "Compile with HDF5" )
	add_executable( MOSGWA ${core_libraries} $<TARGET_OBJECTS:io_hdf5> )
	target_link_libraries( MOSGWA ${external_core_libraries} ${BLAS_LIBRARIES} ${HDF5_LIBRARIES} )
else ( HDF5_FOUND )
	message( STATUS "Compile without HDF5" )
	add_executable( MOSGWA ${core_libraries} )
	target_link_libraries( MOSGWA ${external_core_libraries} ${BLAS_LIBRARIES} )
endif ( HDF5_FOUND )

add_custom_command(
  TARGET MOSGWA
  POST_BUILD
  COMMAND cmake -E remove timestamp.txt
  VERBATIM
)

add_library( TestSuite OBJECT TestSuite.cpp TestSuite.hpp )
add_dependencies( TestSuite TestSuite.tpl )
foreach ( unittest ExceptionTest FirthimizerTest QRuncherTest )
	add_executable( ${unittest} ${unittest}.cpp $<TARGET_OBJECTS:TestSuite> $<TARGET_OBJECTS:main> $<TARGET_OBJECTS:io> $<TARGET_OBJECTS:io_hdf5> $<TARGET_OBJECTS:linalg> $<TARGET_OBJECTS:search> $<TARGET_OBJECTS:imputation> $<TARGET_OBJECTS:parser> $<TARGET_OBJECTS:valuation> $<TARGET_OBJECTS:util> $<TARGET_OBJECTS:minimization> $<TARGET_OBJECTS:lookup> )
	target_link_libraries( ${unittest} unit++ blas gsl hdf5 gomp z pthread )
	add_test( ${unittest} ${unittest} )
endforeach ( unittest )
