#---------------------------------------------------------------------------------
# TARGET is the name of the output
# BUILD is the directory where object files & intermediate files will be placed
# SOURCES is a list of directories containing source code
# INCLUDES is a list of directories containing extra header files
#---------------------------------------------------------------------------------


## top directory where everything happens
TOPDIR   := $(shell pwd)
TARGET   := $(TOPDIR)/bin
BUILD    := $(TOPDIR)/build
SRC      := $(TOPDIR)/src

CC = gcc
CXX = g++
CFLAGS = -g -O3
CXXFLAGS = $(CFLAGS)

#---------------------------------------------------------------------------------
# any extra libraries we wish to link with the projects
#---------------------------------------------------------------------------------
LIBS  := -lm -lgsl -lgslcblas 


## create directory structure
makeDirs:
	@mkdir -p $(TARGET) $(BUILD)

clean:
	rm -rf $(BUILD)
	rm $(SRC)/*/*.gch

all: makeDirs \
        saint-reformat \
	saint-spc-noctrl-matrix \
	saint-spc-noctrl \
	saint-spc-ctrl \
        saint-int-ctrl \
	clean
	@echo -e "\n\n### All executables are in $(TOPDIR)/bin ###\n\n"

saint-reformat:
	$(CC) $(CFLAGS) -c $(SRC)/SAINTreformat/*.c $(SRC)/SAINTreformat/*.h
	mv *.o $(BUILD)
	$(CC) $(LIBS) $(BUILD)/*.o -o $(TARGET)/saint-reformat
	rm $(BUILD)/*.o
	@echo
	@echo

	
saint-spc-noctrl-matrix:
	$(CC) $(CFLAGS) -c $(SRC)/SAINTspc-noctrl-matrix/*.c $(SRC)/SAINTspc-noctrl-matrix/*.h
	mv *.o $(BUILD)
	$(CC) $(LIBS) $(BUILD)/*.o -o $(TARGET)/saint-spc-noctrl-matrix
	rm $(BUILD)/*.o
	@echo
	@echo

saint-spc-noctrl:
	$(CC) $(CFLAGS) -c $(SRC)/SAINTspc-noctrl/*.c $(SRC)/SAINTspc-noctrl/*.h
	mv *.o $(BUILD)
	$(CC) $(LIBS) $(BUILD)/*.o -o $(TARGET)/saint-spc-noctrl
	rm $(BUILD)/*.o
	@echo
	@echo


saint-spc-ctrl:
	$(CC) $(CFLAGS) -c $(SRC)/SAINTspc-ctrl/*.c $(SRC)/SAINTspc-ctrl/*.h
	mv *.o $(BUILD)
	$(CC) $(LIBS) $(BUILD)/*.o -o $(TARGET)/saint-spc-ctrl
	rm $(BUILD)/*.o
	@echo
	@echo


saint-int-ctrl:
	$(CC) $(CFLAGS) -c $(SRC)/SAINTint-ctrl/*.c $(SRC)/SAINTint-ctrl/*.h
	mv *.o $(BUILD)
	$(CC) $(LIBS) $(BUILD)/*.o -o $(TARGET)/saint-int-ctrl
	rm $(BUILD)/*.o
	@echo
	@echo




