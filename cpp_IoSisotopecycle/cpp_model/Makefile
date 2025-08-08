CXX = g++
CXXFLAGS = -std=c++17 -O2
INCLUDES = -Iheaders

all: run_model

run_model: running_model.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o run_model running_model.cpp

clean:
	rm -f run_model