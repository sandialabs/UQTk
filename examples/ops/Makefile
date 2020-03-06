SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET  = Ops.x
TARGET2 = LogComp.x

SRCS  = Ops.cpp Utils.cpp
SRCS2 = LogComp.cpp Utils.cpp

OBJS  = $(SRCS:.cpp=.o)
OBJS2 = $(SRCS2:.cpp=.o)

INCS = Utils.h
LIBINCDIR = $(TOPDIR)/src_cpp/include
LIBBINDIR = $(TOPDIR)/src_cpp/lib
INCDIRS = -I. -I$(LIBINCDIR)
LIBS = -L$(LIBBINDIR) -lUQTk -lquad -luqtktools -lcvode-2.6.0 -ldsfmt -l$(LAPACK) -l$(SLATEC) -l$(BLAS) $(FCLIB)
LDEP = $(LIBBINDIR)/libUQTk.a $(LIBBINDIR)/libquad.a $(LIBBINDIR)/libuqtktools.a $(LIBBINDIR)/libcvode-2.6.0.a \
       $(LIBBINDIR)/libdsfmt.a $(LIBBINDIR)/lib$(LAPACK).a $(LIBBINDIR)/lib$(SLATEC).a $(LIBBINDIR)/lib$(BLAS).a
DEFS = -D__$(FTNNAME)
.SUFFIXES : .cpp

all: $(TARGET) $(TARGET2)

$(TARGET): $(OBJS) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET) $(OBJS) $(LIBS)

$(TARGET2): $(OBJS2) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET2) $(OBJS2) $(LIBS)

clean:
	rm -f $(OBJS) $(OBJS2) $(TARGET) $(TARGET2)
	rm -f samples.*.dat samples.*.dat.pdf

.f.o:
	$(F77) $(FFLAGS) $(DEFS) $(INCDIRS) -c $*.f

.c.o:
	$(CC) $(CCFLAGS) $(DEFS) $(INCDIRS) -c $*.c

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

examples: $(TARGET)
	@echo "Examples of Elementary PC operations:"
	./$(TARGET)
	./plot_pdf.py samples.a.dat
	./plot_pdf.py samples.loga.dat

test: $(TARGET)
	@echo "The test target has not been implemented for $(TARGET)"
