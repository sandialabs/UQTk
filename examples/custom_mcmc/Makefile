UQTK_INS = ../../../UQTk-install
FORTRAN_LIBLOC = /opt/local/lib/gcc5
CXX = /opt/local/bin/g++-mp-5
CXXFLAGS = -O2 
#-std=c++11
TARGET = custom_inf
#TARGET = custom_sam

SRCS   = $(TARGET).cpp
OBJS   = $(SRCS:.cpp=.o)
INCS   = 

LIBINCDIR1 = $(UQTK_INS)/include/uqtk
LIBINCDIR2 = $(UQTK_INS)/include/dep
INCDIRS    = -I. -I$(LIBINCDIR1) -I$(LIBINCDIR2)

LIBDIR = $(UQTK_INS)/lib
LIBSU  = -L$(LIBDIR) -luqtk 
LIBSD  = -L$(LIBDIR) -ldepcvode -ldepnvec -ldepuqtk
FLIB   = -L$(FORTRAN_LIBLOC) -lgfortran
LIBS   =  $(LIBSU) $(LIBSD) $(FLIB) 
LDEP   = $(LIBDIR)/libuqtk.a 

DEFS = -D__wsu $(DFLAGS)




.SUFFIXES : .cpp

all: $(TARGET)

$(TARGET): $(OBJS) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET) *.dat 

.cpp.o: $(INCS) 
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

