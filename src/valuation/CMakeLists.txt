add_library( valuation OBJECT CachingValuation.cpp CachingValuation.hpp Valuation.cpp Valuation.hpp )

add_library( TestValuation OBJECT TestValuation.cpp TestValuation.hpp )
foreach ( unittest CachingValuationTest )
	add_executable( ${unittest} ${unittest}.cpp $<TARGET_OBJECTS:TestSuite> $<TARGET_OBJECTS:TestValuation> $<TARGET_OBJECTS:valuation> $<TARGET_OBJECTS:util> $<TARGET_OBJECTS:lookup> )
	target_link_libraries( ${unittest} unit++ )
	add_test( ${unittest} ${unittest} )
endforeach ( unittest )
