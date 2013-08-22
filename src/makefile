PARSER_OBS = parser/Parslet.o parser/ParsletTemplate.o parser/NameParslet.o parser/SectionParslet.o parser/VariableParslet.o parser/BoolParslet.o parser/IntParslet.o parser/DoubleParslet.o parser/StringParslet.o parser/VectorStringParslet.o parser/ChoiceParslet.o parser/ConfigParser.o
LINALG_OBS = linalg/package.o linalg/Vector.o linalg/AutoVector.o linalg/Matrix.o linalg/AutoMatrix.o linalg/Permutation.o linalg/AutoPermutation.o
IO_OBS = io/Input.o io/InputAdapter.o io/PlinkInput.o io/Hdf5Input.o io/Hdf5Constants.o
HDF5_OBS = io/hdf5/Id.o io/hdf5/File.o io/hdf5/Datatype.o io/hdf5/Dataset.o io/hdf5/StringList.o io/hdf5/DoubleList.o io/hdf5/DoubleTable.o
MINIMIZATION_OBS = minimization/GslMinimizerHolder.o minimization/GslMinimizerRetriever.o minimization/Minimizable.o minimization/Minimizer.o
LOOKUP_OBS = lookup/package.o lookup/ModelIndex.o lookup/MutableModelIndex.o
OBJ = $(PARSER_OBS) $(LINALG_OBS) $(IO_OBS) $(HDF5_OBS) $(MINIMIZATION_OBS) $(LOOKUP_OBS)
OBJSL =Helpfull.o SortVec.o Parameter.o Log.o SNP.o Individual.o MData.o Model.o GenotypeFreq.o QRuncher.o Firthimizer.o ScoreTestShortcut.o GA.o Pool.o Exception.o $(OBJ)

MOBJ = main.o $(OBJSL)
SLIB = Log.hpp Parameter.hpp Helpfull.hpp
ifeq (${HOST}, ims063-067.msi.meduniwien.ac.at) 
     BLAS = -lopenblas -lpthread	
else
     BLAS = -lgslcblas
endif

CC = g++
DEBUG = -g -pg
CFLAGS = -c -Wall -O $(DEBUG)
LFLAGS = -Wall $(DEBUG)
date = $(shell date "+%x")

all: MOSGWA

mosgwa_src.zip:	*.*pp */*.*pp *.h *.tpl makefile
	printf 'string TIP("MOSGWA/HDF5 %s");' `date +%Y-%m-%dT%H-%M-%S` > TIP
	zip -DJX $@ *.*pp */*.*pp *.h *.tpl makefile TIP -x search/* imputation/* Test* *Test.*pp */Test* */*Test.*pp *test.*pp firthAlgorithmCompare.cpp main.cpp main_artur_problem.cpp main_fregene.cpp main_ga.cpp main_hubbart.cpp GslMinimizerAutopointer.*pp

MOSGWA : $(MOBJ)
	$(CC) -static $(LFLAGS) $(MOBJ) -o MOSGWA -lgsl $(BLAS) -lgomp -lhdf5 -lz -lpthread -lm

%.o :	%.cpp %.hpp
	$(CC) $(CFLAGS) -o $@ $<

main.o : main.cpp  MData.cpp MData.hpp Model.cpp Model.hpp $(SLIB)
	$(CC) $(CFLAGS) main.cpp

Parameter.o : Parameter.hpp Parameter.cpp parser/ConfigParser.hpp 
	$(CC) $(CFLAGS) Parameter.cpp

MData.o : MData.cpp MData.hpp Individual.hpp GenotypeFreq.hpp PermSort.hpp  $(SLIB) 
	$(CC) $(CFLAGS) MData.cpp

Model.o : Model.cpp Model.hpp GenotypeFreq.hpp QRuncher.hpp $(SLIB)
	$(CC) $(CFLAGS) Model.cpp

GenotypeFreq.o : GenotypeFreq.cpp GenotypeFreq.hpp MData.hpp
	$(CC) $(CFLAGS) GenotypeFreq.cpp

clean:
	rm *.o */*.o */*/*.o  MOSGWA

