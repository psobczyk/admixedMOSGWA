add_library( io_hdf5 OBJECT
	Dataset.hpp
	Dataset.cpp
	Dataspace.hpp
	Datatype.hpp
	Datatype.cpp
	DoubleList.hpp
	DoubleList.cpp
	DoubleTable.hpp
	DoubleTable.cpp
	File.hpp
	File.cpp
	Id.hpp
	Id.cpp
	Object.hpp
	StringList.hpp
	StringList.cpp
)

add_dependencies( io_hdf5
	Dataspace.tpl
	Object.tpl
)
