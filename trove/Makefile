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
INCLUDES	= -I../src -I../trove
VPATH		= ../src:../test
LIB		= -lpthread

CXXFLAGS	= $(OPTS) $(LANG) $(PICKY) $(INCLUDES)

TESTS		= testmult.cpp testvec.cpp testcsr.cpp testcoo.cpp testsparse.cpp bench.cpp
SOURCES		= AOSMatrix.cpp COOMatrix.cpp CSRMatrix.cpp Matrix.cpp Vector.cpp

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
	@ $(CXX) -MM $(LANG) $(INCLUDES) $(SOURCES) > makedep
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
AOSMatrix.o: AOSMatrix.cpp AOSMatrix.hpp Vector.hpp
COOMatrix.o: COOMatrix.cpp COOMatrix.hpp Vector.hpp
CSRMatrix.o: CSRMatrix.cpp CSRMatrix.hpp Vector.hpp
Matrix.o: Matrix.cpp Matrix.hpp Vector.hpp
Vector.o: Vector.cpp Vector.hpp
# DEPENDENCIES MUST END AT END OF FILE
# IF YOU PUT STUFF HERE IT WILL GO AWAY
# see make depend above
