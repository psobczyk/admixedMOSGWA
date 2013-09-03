add_library( util OBJECT Array.hpp Cache.cpp Cache.hpp ReferenceCounted.cpp ReferenceCounted.hpp )
add_dependencies( util ${util_deps} Array.tpl )

foreach ( unittest ReferenceCountedTest ArrayTest CacheTest )
	add_executable( ${unittest} ${unittest}.cpp $<TARGET_OBJECTS:TestSuite> $<TARGET_OBJECTS:util> )
	target_link_libraries( ${unittest} unit++ )
	add_test( ${unittest} ${unittest} )
endforeach ( unittest )
