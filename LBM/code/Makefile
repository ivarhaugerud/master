CXXFLAGS = -Wall -O3 -larmadillo -std=c++11

program: MainClass.o test_heat_both_ways.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm program *.o
