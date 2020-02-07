SHELL = /bin/bash
TARGET = dfi_line.exe

UQTK_LIBS= -L$(UQTK_INS)/lib \
-luqtk -luqtkdfi -ldepuqtk -l sundials_cvode -lblas -llapack -lgfortran
UQTK_INCL= -I$(UQTK_INS)/include/uqtk -I$(UQTK_INS)/include/dep -I$(UQTK_INS)/include

CXX = g++
CFLAGS= -std=c++11
all: $(TARGET)

userFunctions_line.o: userFunctions_line.cpp
	$(info compile the user functions...)
	$(CXX) $(CFLAGS) $(UQTK_INCL) -c userFunctions_line.cpp $(UQTK_LIBS)
	ar -crs libUserFunc_line.a userFunctions_line.o	

main.o: main.cpp
	$(info make the driver...)
	$(CXX) $(CFLAGS) $(UQTK_INCL) -c main.cpp $(UQTK_LIBS)

$(TARGET): userFunctions_line.o main.o
	$(info build the dfi executeable...)
	$(CXX) $(UQTK_INCL) $(CFLAGS) -o $(TARGET) main.o $(UQTK_LIBS) libUserFunc_line.a  
clean:
	rm -f *.a *.o $(TARGET) 
	rm -rf html latex 
