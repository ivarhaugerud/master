CXXFLAGS = -Wall -O3 -larmadillo -std=c++11

program: MainClass.o injection_peaks.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm program *.o
