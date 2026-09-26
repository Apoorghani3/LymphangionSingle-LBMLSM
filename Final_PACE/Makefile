# LBLS linux Makefile

SHELL = /bin/sh
#.IGNORE:

# Definitions

ROOT =	lbls
EXE =	$(ROOT)
SRC =	$(wildcard source/*.cpp)
INC =	$(wildcard source/*.h)
OBJ = 	$(SRC:.cpp=.o)

# System-specific settings

CC = 		g++
CCFLAGS =	-w -fpermissive -O3 -fopenmp
DEPFLAGS =	-M
LINK =		g++
LINKFLAGS =	-O3 -fopenmp
USRLIB =	
SYSLIB =	
SIZE =		size


.DEFAULT:

# Link rule

$(EXE):	$(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(USRLIB) $(SYSLIB) -o $(EXE)
	$(SIZE) $(EXE)

# Compilation rules

%.o:%.cpp $(INC)
	$(CC) $(CCFLAGS) -c $< -o $@

#%.d:%.cpp
#	$(CC) $(CCFLAGS) $(DEPFLAGS) $< > $@

# Individual dependencies

#DEPENDS = $(OBJ:.o=.d)
#include $(DEPENDS)


#makefile.dep:$(SRC) $(INC)
#	for i in $(SRC); do $(CC) $(CCFLAGS) $(DEPFLAGS) "$${i}"; done > $@

#include makefile.dep




