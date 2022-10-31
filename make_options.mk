# User Configurable Options
PROJECT_ROOT=/Users/amedhi/Projects/Codes/VMC/MainVMC

#-------------------------------------------------------------
# need mpi version
#MPI=HAVE_BOOST_MPI

#-------------------------------------------------------------
# 1. Set compiler option
CXX=clang++ -std=c++11 # Clang compiler 
#CXX=g++ -std=c++11 # GNU GCC compiler

ifeq ($(MPI), HAVE_BOOST_MPI)
  CXX=/opt/openmpi/bin/mpicxx -std=c++11
  CPPFLAGS=-D$(MPI)
else
  CPPFLAGS=
endif

#-------------------------------------------------------------
WAVEFUNC=REAL
ifeq ($(WAVEFUNC), REAL)
  CPPFLAGS += -DREAL_WAVEFUNCTION -DEIGEN_NO_DEBUG 
else
  CPPFLAGS += -DEIGEN_NO_DEBUG 
endif
#EIGEN_USE_MKL=EIGEN_USE_MKL_ALL
ifeq ($(EIGEN_USE_MKL), USE_INTEL_MKL_ALL)
  CPPFLAGS += -DEIGEN_USE_MKL_ALL
endif


#-------------------------------------------------------------
# 2. Compile flags 
# Flags to give the compiler for "release mode"
OPTFLAGS= -Wall -O3
#OPTFLAGS=-Wall -pg
#OPTFLAGS=-DDEBUG_MODE -g -Wall -pedantic
# Flags to give the compiler for "debug mode"
#DEBUGFLAGS=-DDEBUG_MODE -g -Wall -pedantic

#-------------------------------------------------------------
# 3. Boost and Eigen library
# Flags to give the compiler for "release mode"
BOOST_INCLUDE=-I/usr/local/include
EIGEN_INCLUDE=-I/usr/local/include/Eigen 

# Boost MPI library
ifeq ($(MPI), HAVE_BOOST_MPI)
  BOOST_LIBS=-lboost_filesystem -lboost_system -lboost_mpi -lboost_serialization
  BOOST_LDFLAGS=-L/usr/local/lib
else 
  BOOST_LIBS=-lboost_filesystem -lboost_system
  BOOST_LDFLAGS=-L/usr/local/lib
endif

MKL_INCLUDE=-I/opt/intel/mkl/include/intel64/lp64
MKL_LDFLAGS=-L/opt/intel/mkl/lib
MKL_LIBS=-lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#NLOPT_INCLUDE=-I/Users/amedhi/projects/Codes/vmc++/libs/include
#NLOPT_LDFLAGS=-L/Users/amedhi/projects/Codes/vmc++/libs/lib
#NLOPT_LIBS=-lnlopt

INCLUDE = $(BOOST_INCLUDE) #$(MKL_INCLUDE)
ifneq ($(BOOST_INCLUDE), $(EIGEN_INCLUDE))
  INCLUDE += $(EIGEN_INCLUDE)
endif

CXXFLAGS = $(CPPFLAGS) $(OPTFLAGS) $(INCLUDE) #$(NLOPT_INCLUDE)
LDFLAGS = $(BOOST_LDFLAGS) $(CERES_LDFLAG) #$(NLOPT_LDFLAGS)  
LIBS = $(BOOST_LIBS) #$(NLOPT_LIBS)  

ifeq ($(EIGEN_USE_MKL), USE_INTEL_MKL_ALL)
  INCLUDE += $(MKL_INCLUDE)
  LDFLAGS += $(MKL_LDFLAGS)
  LIBS += $(MKL_LIBS) 
endif

#-------------------------------------------------------------
# 4. Build direcory
PREFIX=$(PROJECT_ROOT)
BUILD_DIR=$(PREFIX)/build

VMC_LIBDIR=$(PREFIX)/lib
VMC_INCLUDE=$(PREFIX)/include
#VMC_CXXFLAGS= $(VMC_OPTFLAGS) $(INCLUDE) -I$(VMC_INCLUDE)
#VMC_LDFLAGS=$(BOOST_LDFLAGS) -L$(VMC_LIBDIR)
#VMC_LIBS=$(BOOST_LIBS) -lvmc++
#-------------------------------------------------------------
