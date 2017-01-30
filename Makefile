#-------------------------------------------------------------
# Makefile for cmc++ library
#-------------------------------------------------------------
#include ../root_dir.mk # must be included first
include ./options.mk
#-------------------------------------------------------------
# Source files
SRCS = scheduler/cmdargs.cc 
SRCS+= scheduler/inputparams.cc 
SRCS+= scheduler/taskparams.cc 
SRCS+= scheduler/worker.cc 
SRCS+= scheduler/master_scheduler.cc
SRCS+= scheduler/scheduler.cc
SRCS+= expression/expression.cc 
SRCS+= expression/tokens.cc 
SRCS+= expression/functions.cc
SRCS+= expression/objects.cc
SRCS+= lattice/lattice.cc
SRCS+= lattice/latticelibrary.cc
SRCS+= lattice/graph.cc
#SRCS+= model/qn.cc
#SRCS+= model/quantum_operator.cc
#SRCS+= model/sitebasis.cc
SRCS+= model/hamiltonian_term.cc
SRCS+= model/model.cc
SRCS+= model/modellibrary.cc
SRCS+= variational/wavefunction.cc
#SRCS+= observable/observables.cc
SRCS+= montecarlo/random.cc
#SRCS+= montecarlo/observable_operator.cc
#SRCS+= montecarlo/measurement.cc
SRCS+= montecarlo/simulator.cc
SRCS+= main.cc
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
         model/modelparams.h  model/hamiltonian_term.h model/model.h \
	 variational/wavefunction.h \
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
OBJS=$(patsubst %.cc,$(BUILD_DIR)/%.o,$(VMC_SRCS))
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

%.o: %.cc
	$(VMC_CXX) -c $(CXXBFLAGS) -o $@ $<


# Include all .d files
-include $(DEPS)

$(BUILD_DIR)/%.o: %.cc
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
