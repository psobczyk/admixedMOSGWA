add_library( io_hdf5 OBJECT Dataset.cpp Dataset.hpp Dataspace.hpp Datatype.cpp Datatype.hpp DoubleList.cpp DoubleList.hpp DoubleTable.cpp DoubleTable.hpp File.cpp File.hpp Id.cpp Id.hpp Object.hpp StringList.cpp StringList.hpp )
add_dependencies( io_hdf5 Dataspace.tpl Object.tpl )
