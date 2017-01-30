# User Configurable Options
PROJECT_ROOT=/Users/amedhi/projects/Codes/vmc++

# Building the cmc++ lib
#-------------------------------------------------------------
# 1. Set compiler option
VMC_CXX=clang++ -std=c++11 # Clang compiler 
#VMC_CXX=g++ -std=c++11 # GNU GCC compiler

ifndef VMC_CXX
$(error Makefile variable VMC_CXX not defined in options.mk; please define it.)
endif
#-------------------------------------------------------------
# 2. Compile flag options
# Flags to give the compiler for "release mode"
VMC_OPTFLAGS=-Wall -O3
# Flags to give the compiler for "debug mode"
VMC_DEBUGFLAGS=-DDEBUG_MODE -g -Wall -pedantic

#-------------------------------------------------------------
# 3. Boost and Eigen library
# Flags to give the compiler for "release mode"
BOOST_INCLUDE=-I/usr/local/include
EIGEN_INCLUDE=-I/usr/local/include
#BOOST_LIBS=-lboost_system -lboost_filesystem
#BOOST_LDFLAGS=-L/usr/local/lib

INCLUDE = $(BOOST_INCLUDE)
ifneq ($(BOOST_INCLUDE), $(EIGEN_INCLUDE))
INCLUDE += $(EIGEN_INCLUDE)
endif

VMC_CXXBFLAGS= $(VMC_OPTFLAGS) $(INCLUDE)
#-------------------------------------------------------------
# 4. Where to put the 'cmc' library & the includes
PREFIX=$(PROJECT_ROOT)
BUILD_DIR=$(PREFIX)/build
VMC_LIBDIR=$(PREFIX)/lib
VMC_INCLUDE=$(PREFIX)/include
VMC_CXXFLAGS= $(VMC_OPTFLAGS) $(INCLUDE) -I$(VMC_INCLUDE)
VMC_LDFLAGS=$(BOOST_LDFLAGS) -L$(VMC_LIBDIR)
VMC_LIBS=$(BOOST_LIBS) -lcmc++
#-------------------------------------------------------------
