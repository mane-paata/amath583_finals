#
# This file is part of the course materials for AMATH483/583 at the University of Washington,
# Spring 2018
#
# Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
# https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Author: Andrew Lumsdaine
# 

LANG		= -std=c++14
PICKY		= -Wall

CXX		= c++
OPTS		= -Ofast -march=native -DNDEBUG 
REPORTING   	= -Rpass=.*
INCLUDES	= -I../src -I../test -I../trove 
VPATH		= ../src:../test:../trove
LIB		= -lpthread

CXXFLAGS	= $(OPTS) $(LANG) $(PICKY) $(INCLUDES) 

TESTS		= testmult.cpp testvec.cpp testcsr.cpp testcoo.cpp testsparse.cpp \
		  bench.cpp sparsebench.cpp ir_driver.cpp jacobi_driver.cpp
SOURCES		= AOSMatrix.cpp COOMatrix.cpp CSRMatrix.cpp Matrix.cpp Vector.cpp \
		  Grid.cpp Stencil.cpp ir.cpp

HEADERS		= $(SOURCES:.cpp=.hpp) 
OBJECTS		= $(SOURCES:.cpp=.o) $(TESTS:.cpp=.o)
TARGETS		= $(TESTS:.cpp=.exe)
PCHS		= $(HEADERS:=.gch)



.PHONY		: defreport optreport clean depend all

all		: $(TARGETS)

%.exe        	: %.o
		  $(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)

%.o 		: %.cpp
		  $(CXX) -c $(CXXFLAGS) $< -o $@

%.s 		: %.cpp
		  $(CXX) -S $(CXXFLAGS) $<

$(TARGETS) 	: $(SOURCES:.cpp=.o)

defreport	:
		  $(CXX) -dM -E -x c++ /dev/null

optreport	:
		  echo 'int;' | $(CXX) -xc++ $(CXXFLAGS) - -o /dev/null -\#\#\#

clean		:
		  /bin/rm -f $(TARGETS) $(OBJECTS) $(PCHS) Matrix.s a.out *~

depend: $(ALL_SRC)
	@ $(CXX) -MM $(LANG) $(INCLUDES) $(TESTS) > makedep
	@ echo '/^# DO NOT DELETE THIS LINE/+2,$$d' >eddep
	@ echo '$$r makedep' >>eddep
	@ echo 'w' >>eddep
	@ cp Makefile Makefile.bak
	@ ed - Makefile < eddep
	@ /bin/rm eddep makedep
	@ echo '# DEPENDENCIES MUST END AT END OF FILE' >> Makefile
	@ echo '# IF YOU PUT STUFF HERE IT WILL GO AWAY' >> Makefile
	@ echo '# see make depend above' >> Makefile

# The following 4 (yes 4) lines must be in all the subdirectory makefiles
#-----------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend uses it
# DEPENDENCIES MUST END AT END OF FILE
testmult.o: testmult.cpp ../trove/Matrix.hpp ../trove/Vector.hpp
testvec.o: testvec.cpp ../trove/Vector.hpp
testcsr.o: testcsr.cpp ../trove/CSRMatrix.hpp ../trove/Vector.hpp \
 ../trove/Vector.hpp
testcoo.o: testcoo.cpp ../trove/COOMatrix.hpp ../trove/Vector.hpp \
 ../trove/Vector.hpp
testsparse.o: testsparse.cpp ../trove/AOSMatrix.hpp ../trove/Vector.hpp \
 ../trove/COOMatrix.hpp ../trove/CSRMatrix.hpp ../trove/Matrix.hpp \
 ../trove/Timer.hpp ../trove/Vector.hpp
bench.o: bench.cpp ../trove/Matrix.hpp ../trove/Vector.hpp \
 ../trove/Timer.hpp
sparsebench.o: sparsebench.cpp ../trove/COOMatrix.hpp ../trove/Vector.hpp \
 ../trove/Matrix.hpp ../trove/Timer.hpp ../trove/Vector.hpp
ir_driver.o: ir_driver.cpp ../trove/Grid.hpp ../trove/Stencil.hpp \
 ../trove/Grid.hpp ../trove/ir.hpp ../trove/Stencil.hpp
jacobi_driver.o: jacobi_driver.cpp ../trove/Grid.hpp
# DEPENDENCIES MUST END AT END OF FILE
# IF YOU PUT STUFF HERE IT WILL GO AWAY
# see make depend above
