CXXFLAGS = -Wall -O3 -larmadillo -std=c++11

program: MainClass.o benchmarking.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm program *.o
