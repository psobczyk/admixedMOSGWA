PARSER_OBS = parser/Parslet.o parser/ParsletTemplate.o parser/NameParslet.o parser/SectionParslet.o parser/VariableParslet.o parser/BoolParslet.o parser/IntParslet.o parser/DoubleParslet.o parser/StringParslet.o parser/VectorStringParslet.o parser/ChoiceParslet.o parser/ConfigParser.o
LINALG_OBS = linalg/package.o linalg/Vector.o linalg/AutoVector.o linalg/Matrix.o linalg/AutoMatrix.o linalg/Permutation.o linalg/AutoPermutation.o
IO_OBS = io/Input.o io/PlinkInput.o io/Hdf5Input.o
MINIMIZATION_OBS = minimization/GslMinimizerHolder.o minimization/GslMinimizerRetriever.o minimization/Minimizable.o minimization/Minimizer.o
LOOKUP_OBS = lookup/package.o lookup/ModelIndex.o lookup/MutableModelIndex.o
OBJ = $(PARSER_OBS) $(LINALG_OBS) $(IO_OBS) $(MINIMIZATION_OBS) $(LOOKUP_OBS)
OBJSL =Helpfull.o SortVec.o Parameter.o Log.o SNP.o Individual.o MData.o Model.o GenotypeFreq.o QRuncher.o Firthimizer.o ScoreTestShortcut.o GA.o Pool.o Exception.o $(OBJ)

#SortVec.o Helpfull.o Parameter.o Log.o SNP.o MData.o  Model.o  GenotypeFreq.o $(OBJ)
MOBJ = main.o $(OBJSL)
FOBJ = main_fregene.o $(OBJSL)
EOBJ = main_erich.o $(OBJSL)
AOBJ = main_artur_problem.o $(OBJSL)
GAOBJ = main_ga.o $(OBJSL)
SLIB = Log.hpp Parameter.hpp Helpfull.hpp
ifeq (${HOST}, ims063-067.msi.meduniwien.ac.at) 
     BLAS = -lopenblas -lpthread	
else
     BLAS = -lgslcblas
endif

CC = g++ -std=c++0x -lhdf5 -lz -I/usr/lib/openmpi/include
CC = h5c++ -echo -std=c++0x -I/usr/lib/openmpi/include
DEBUG = -g -pg
CFLAGS =  -Wall    -Wno-sign-compare -c -fopenmp -O $(DEBUG)
CFLAGS =  -static -Wall    -Wno-sign-compare -c -O $(DEBUG)
#CFLAGS = -static -Wall  -Wno-sign-compare -c -I/usr/local/lib/ -O $(DEBUG)
LFLAGS = -Wall -L/usr/local/lib64/ -fopenmp $(DEBUG)
LFLAGS = -Wall -L/usr/local/lib64/ $(DEBUG)
ifeq (${HOST},fs) #utmost oldfashioned
     #CC = h5c++ 
     CFLAGS =  -Wall  -Wno-deprecated-declarations  -Wno-sign-compare -c  -I/home/ffrommlet/erich/include  -O $(DEBUG)
     LFLAGS = -Wall -L/usr/local/lib64/  -L/home/ffrommlet/erich/lib $(DEBUG) 
endif
date = $(shell date "+%x")
#-lgfortran 

all: MOSGWA

mosgwa_src.zip:	*.*pp */*.*pp *.h *.tpl makefile
	printf 'string TIP("MOSGWA/HDF5 %s");' `date +%Y-%m-%dT%H-%M-%S` > TIP
	zip -DJX $@ *.*pp */*.*pp *.h *.tpl makefile TIP -x search/* imputation/* Test* *Test.*pp */Test* */*Test.*pp *test.*pp firthAlgorithmCompare.cpp main.cpp main_artur_problem.cpp main_fregene.cpp main_ga.cpp main_hubbart.cpp GslMinimizerAutopointer.*pp

MOSGWA : $(MOBJ)
	$(CC) -static $(LFLAGS) $(MOBJ) -o MOSGWA -lgsl $(BLAS) -lpthread -lm -lhdf5 -lgomp

MOSGWA_erich : $(EOBJ)
	#$(CC)   $(LFLAGS) $(EOBJ) -o MOSGWA -lgsl -lgslcblas  -lm
	$(CC) -static   $(LFLAGS) $(EOBJ)   -o MOSGWA -lgsl $(BLAS) -lpthread -lm -lhdf5 -lgomp

MOSGWA_fregene : $(FOBJ)
	$(CC)   $(LFLAGS) $(FOBJ) -o MOSGWA_fregene  -lgsl -lgslcblas -lm #-lgfortran
A:$(AOBJ)
	$(CC)   $(LFLAGS) $(AOBJ)   -o A -lgsl $(BLAS) -lpthread -lm -lhdf5 -lgfortran  

MOSGWA_GA : $(GAOBJ)
	$(CC) $(LFLAGS) $(GAOBJ) -o GA_MOSGWA -lgsl -lgslcblas -lm

%.o :	%.cpp %.hpp
	$(CC) $(CFLAGS) -o $@ $<

main.o : main.cpp  MData.cpp MData.hpp Model.cpp Model.hpp $(SLIB)
	$(CC) $(CFLAGS) main.cpp

main_ga.o : main_ga.cpp MData.hpp Model.hpp GA.hpp

main_erich.o : main_erich.cpp  MData.cpp MData.hpp  Model.cpp Model.hpp   $(SLIB)
	$(CC) $(CFLAGS) main_erich.cpp

main_fregene.o : main_fregene.cpp MData.hpp  Model.hpp $(SLIB)
	$(CC) $(CFLAGS)  main_fregene.cpp 

Parameter.o : Parameter.hpp Parameter.cpp parser/ConfigParser.hpp 
	$(CC) $(CFLAGS) Parameter.cpp
main_artur_problem.o: main_artur_problem.cpp MData.hpp  Model.hpp $(SLIB)
	$(CC) $(CFLAGS) main_artur_problem.cpp
MData.o : MData.cpp MData.hpp OneSnpAllInd.hpp Individual.hpp GenotypeFreq.hpp PermSort.hpp  $(SLIB) 
	$(CC) $(CFLAGS) MData.cpp

Model.o : Model.cpp Model.hpp GenotypeFreq.hpp QRuncher.hpp $(SLIB)
	$(CC) $(CFLAGS) Model.cpp

GenotypeFreq.o : GenotypeFreq.cpp GenotypeFreq.hpp MData.hpp
	$(CC) $(CFLAGS) GenotypeFreq.cpp

clean:
	rm *.o */*.o  MOSGWA

