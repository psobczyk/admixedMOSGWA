add_library( util OBJECT
	Array.hpp
	Cache.hpp
	Cache.cpp
	ReferenceCounted.hpp
	ReferenceCounted.cpp
	LogFactorial.hpp
	LogFactorial.cpp
)

add_dependencies( util
	Array.tpl
)
