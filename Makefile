#-------------------------------------------------------------
# Makefile for cmc++ library
#-------------------------------------------------------------
#include ../root_dir.mk # must be included first
include ./options.mk
#-------------------------------------------------------------
# Source files
SRCS = scheduler/cmdargs.cpp 
SRCS+= scheduler/inputparams.cpp 
SRCS+= scheduler/taskparams.cpp 
SRCS+= scheduler/worker.cpp 
SRCS+= scheduler/master_scheduler.cpp
SRCS+= scheduler/scheduler.cpp
#SRCS+= xml/pugixml.cpp 
SRCS+= expression/expression.cpp 
SRCS+= expression/tokens.cpp 
SRCS+= expression/functions.cpp
SRCS+= expression/objects.cpp
SRCS+= lattice/lattice.cpp
SRCS+= lattice/latticelibrary.cpp
SRCS+= lattice/graph.cpp
#SRCS+= model/qn.cpp
#SRCS+= model/quantum_operator.cpp
#SRCS+= model/sitebasis.cpp
SRCS+= model/hamiltonian_term.cpp
SRCS+= model/model.cpp
SRCS+= model/modellibrary.cpp
SRCS+= variational/blochbasis.cpp
SRCS+= variational/mf_model.cpp
SRCS+= variational/bcsfunction.cpp
SRCS+= variational/wavefunction.cpp
#SRCS+= observable/observables.cpp
SRCS+= montecarlo/random.cpp
#SRCS+= montecarlo/observable_operator.cpp
#SRCS+= montecarlo/measurement.cpp
SRCS+= montecarlo/simulator.cpp
SRCS+= main.cpp
VMC_SRCS = $(addprefix src/,$(SRCS))
#-------------------------------------------------------------
# Headers
HDRS=scheduler/optionparser.h scheduler/cmdargs.h \
         scheduler/inputparams.h scheduler/worker.h scheduler/task.h \
         scheduler/scheduler.h \
         expression/expression.h expression/shunting_yard.h \
         expression/tokens.h expression/functions.h expression/objects.h \
         expression/pack.h \
         lattice/constants.h lattice/lattice.h lattice/graph.h \
	 montecarlo/simulator.h \
         model/modelparams.h  model/quantum_op.h \
	 model/hamiltonian_term.h model/model.h \
	 variational/blochbasis.h \
	 variational/mf_model.h variational/wavefunction.h \
         montecarlo/random.h montecarlo/sitebasisstate.h \
	 simulation.h \
#         model/qn.h model/quantum_operator.h model/sitebasis.h model/modelparams.h \
#         observable/mcdata.h observable/observables.h \
#	 montecarlo/observable_operator.h montecarlo/simulator.h \
	 simulation.h
VMC_HDRS = $(addprefix src/,$(HDRS))
#-------------------------------------------------------------
# Target
TAGT=a.out

# Put all auto generated stuff to this build dir.
ifeq ($(BUILD_DIR), $(CURDIR))
  $(error In-source build is not allowed, choose another build directory)
endif

# All .o files go to BULD_DIR
OBJS=$(patsubst %.cpp,$(BUILD_DIR)/%.o,$(VMC_SRCS))
# GCC/Clang will create these .d files containing dependencies.
DEPS=$(patsubst %.o,%.d,$(OBJS)) 
# compiler flags
CXXFLAGS = $(VMC_INCLUDEFLAGS) $(OPTIMIZATIONS) 

#VMC_INCLDIR=$(VMC_INCLUDE)/cmc
#INCL_HDRS = $(addprefix $(VMC_INCLDIR)/, $(VMC_HDRS))

.PHONY: all
all: $(TAGT) #$(INCL_HDRS)

$(TAGT): $(OBJS)
	$(VMC_CXX) -o $(TAGT) $(OBJS) $(BOOST_LDFLAGS) $(BOOST_LIBS) 

%.o: %.cpp
	$(VMC_CXX) -c $(CXXBFLAGS) -o $@ $<


# Include all .d files
-include $(DEPS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	@echo "$(VMC_CXX) -c $(VMC_CXXBFLAGS) -o $(@F) $(<F)"
	@$(VMC_CXX) -MMD -c $(VMC_CXXBFLAGS) -o $@ $<

$(VMC_INCLDIR)/%.h: %.h 
	@mkdir -p $(@D)
	@echo "Copying $< to 'include'" 
	@cp -f $< $@

# installation
#prefix = ../install#/usr/local
#libdir = $(prefix)/lib
#includedir = $(prefix)/include/cmc++

.PHONY: install
install:	
	@echo "Already installed in $(VMC_LIBDIR) and $(VMC_INCLDIR)" 

.PHONY: clean
clean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
	@echo "Removing $(TAGT)"
	@rm -f $(TAGT) 

.PHONY: bclean
bclean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
