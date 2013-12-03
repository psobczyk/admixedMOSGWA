add_library( util OBJECT Array.hpp Cache.cpp Cache.hpp ReferenceCounted.cpp ReferenceCounted.hpp )
add_dependencies( util Array.tpl )
