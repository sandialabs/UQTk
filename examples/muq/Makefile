TARGET = muqPCE2uqtk
SRCS   = $(TARGET).cpp 
OBJS   = $(SRCS:.cpp=.o)

CXX=/usr/local/bin/g++-4.8

CXXFLAGS=-std=c++11 -Wall           \
         -Wno-maybe-uninitialized   \
         -Wno-sign-compare          \
         -Wno-unknown-pragmas       \
         -Wno-unused-variable       \
         -Wno-unused-local-typedefs

#-------------------------------------------------------------------------------------
#  UQTk
#-------------------------------------------------------------------------------------

INCUQTk=-I$(UQTKdir)/include/uqtk -I$(UQTKdir)/include/dep
LIBSUQTk=-L$(UQTKdir)/lib -luqtk -ldepcvode -ldepnvec -ldepuqtk -lgfortran

#-------------------------------------------------------------------------------------
#  MUQ
#-------------------------------------------------------------------------------------
PATH2MUQ=$(HOME)/Projects/muq-install

INCMUQ=-I$(PATH2MUQ)/include                     \
       -I$(PATH2MUQ)/muq_external/include/eigen3 \
       -I$(PATH2MUQ)/muq_external/include

LIBSMUQ=-L$(PATH2MUQ)/lib -L$(PATH2MUQ)/muq_external/lib \
        -lmuqModelling -lmuqPolynomialChaos -lmuqUtilities

LIBSMUQDEP=-L$(PATH2MUQ)/muq_external/lib -lflann_cpp         \
           -lsundials_cvodes -lsundials_idas                  \
           -lsundials_kinsol -lsundials_nvecserial            \
           -lhdf5_hl -lhdf5 -lboost_filesystem -lboost_system \
           -lboost_graph -lboost_serialization -lboost_regex  \
           -lboost_wserialization


LDFLAGS=-Wl,-rpath,$(PATH2MUQ)/lib -Wl,-rpath,$(PATH2MUQ)/muq_external/lib 

DFLAGS=

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBSUQTk) $(LIBSMUQ) $(LIBSMUQDEP) $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(INCUQTk) $(INCMUQ) -c $*.cpp

clean:
	/bin/rm -rf $(TARGET) $(OBJS)
