
#compiler setup
CXX = g++
MPICXX = mpic++
CXXFLAGS = -std=c++14 -O3 -pthread $(MACRO)

COMMON= core/utils.h core/cxxopts.h core/get_time.h
SERIAL= generate linear_eqn_solver_serial
PARALLEL= linear_eqn_solver_parallel linear_eqn_solver_mpi
ALL= $(SERIAL) $(PARALLEL)

all : $(ALL)

$(SERIAL): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(PARALLEL): %: %.cpp
	$(MPICXX) $(CXXFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)
