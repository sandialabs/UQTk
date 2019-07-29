SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site

TARGET1 = SurfRxnNISP.x
TARGET2 = SurfRxnNISP_MC.x
TARGET3 = SurfRxnISP.x
TARGET4 = SurfRxnDet.x
TARGET5 = SurfRxnInfer.x 


SRCS1 = SurfRxnNISP.cpp Utils.cpp
SRCS2 = SurfRxnNISP_MC.cpp Utils.cpp
SRCS3 = SurfRxnISP.cpp Utils.cpp
SRCS4 = SurfRxnDet.cpp
#SRCS5 = SurfRxnInfer.cpp

OBJS1 = $(SRCS1:.cpp=.o)
OBJS2 = $(SRCS2:.cpp=.o)
OBJS3 = $(SRCS3:.cpp=.o)
OBJS4 = $(SRCS4:.cpp=.o)
#OBJS5 = $(SRCS5:.cpp=.o)

INCS = model.h XMLreader.h posterior.h Utils.h

LIBINCDIR = $(TOPDIR)/src_cpp/include
LIBBINDIR = $(TOPDIR)/src_cpp/lib
INCDIRS   = -I. -I$(LIBINCDIR)

LIBS = -L$(LIBBINDIR) -lUQTk -lquad -luqtkmcmc -luqtktools -llbfgs -lcvode-2.6.0 -ldsfmt \
       -l$(LAPACK) -l$(SLATEC) -l$(BLAS) -lxmlutils -lexpat $(FCLIB)
LDEP = $(LIBBINDIR)/libUQTk.a $(LIBBINDIR)/libquad.a $(LIBBINDIR)/liblbfgs.a $(LIBBINDIR)/libbcs.a \
       $(LIBBINDIR)/libuqtktools.a $(LIBBINDIR)/libcvode-2.6.0.a \
       $(LIBBINDIR)/libdsfmt.a $(LIBBINDIR)/lib$(LAPACK).a \
       $(LIBBINDIR)/lib$(SLATEC).a $(LIBBINDIR)/lib$(BLAS).a

DEFS = -D__$(FTNNAME)

.SUFFIXES : .cpp

all: $(TARGET1) $(TARGET2) $(TARGET3)  $(TARGET4) #$(TARGET5)

$(TARGET1): $(OBJS1) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET1) $(OBJS1) $(LIBS)

$(TARGET2): $(OBJS2) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET2) $(OBJS2) $(LIBS)

$(TARGET3): $(OBJS3) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET3) $(OBJS3) $(LIBS)

$(TARGET4): $(OBJS4) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET4) $(OBJS4) $(LIBS)

#$(TARGET5): $(OBJS5) $(INCS) $(LDEP)
#	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET5) $(OBJS5) $(LIBS)

clean:
	rm -f $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2) $(OBJS3) $(TARGET3)
	rm -f $(OBJS4) $(TARGET4) $(OBJS5) $(TARGET5) *~
	rm -f *.dat *.pyc
	rm -f *.parsed.xml
	rm -f prob2_*.pdf prob2_*.eps prob3_*.pdf surf_rxn_*ISP_*.pdf

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

test: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) #$(TARGET5)
	@echo "The test target has not been implemented for $(TARGET)"

examples: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) #$(TARGET5)
	echo "Running problem 2"
	./prob2.py HG 3 ISP NISP NISP_MC
	echo "Running problem 3"
	./prob3.py
#	echo "Running problem 4"
#	./prob4.py
