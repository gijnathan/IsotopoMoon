CXX = g++
CXXFLAGS = -std=c++17 -O2


all: run_model run_model_input

run_model: running_model.cpp headers/isotope_evolution.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o run_model running_model.cpp

run_model_input: running_model_WITH_INPUTFILE.cpp headers/isotope_evolution.hpp headers/config.hpp headers/config_io.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o run_model_input running_model_WITH_INPUTFILE.cpp

clean:
	rm -f run_model run_model_input
