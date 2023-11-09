# Define the C++ compiler
CXX = g++
# Define any compile-time flags
CXXFLAGS = -Wall -std=c++11

# Define the target executable name
TARGET = tsp_nearest_insertion

# All target (default)
all: $(TARGET)

# Link the target executable:
$(TARGET): main.o
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.o

# Compile the source files into object files:
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

# Clean target
clean:
	rm -f *.o $(TARGET)
