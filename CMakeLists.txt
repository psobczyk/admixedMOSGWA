project( MOSGWA C CXX )
cmake_minimum_required( VERSION 2.8.8 )
set( cmake_verbose_makefile on )

add_definitions( -g -pg -Wall )
add_subdirectory( src build )

enable_testing()
