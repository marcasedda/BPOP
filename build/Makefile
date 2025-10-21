#
# Makefile for BPOP
# Manuel Arca Sedda 2025
#

CXX = g++

RESULT   = BPOP.x
CXXFLAGS = -O3 -I../include
LDFLAGS = -o $(RESULT)

VPATH= ../:../src:../include

INC = input_params.h funcGR.h config.h
SOURCE = Main_Bpop_dev.cpp funcGR.cpp config.cpp

.SUFFIXES : .o .cpp

OBJECTS = $(SOURCE:.cpp=.o)

TARGET: $(RESULT) 

$(RESULT): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS)  $(OBJECTS)

$(OBJECTS): $(INC)
