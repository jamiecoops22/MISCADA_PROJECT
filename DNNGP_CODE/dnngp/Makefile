CXX = g++
OPTFLAGS = -O3 -fopenmp
CXXFLAGS = -I/opt/OpenBLAS/include 
LIBS = -I/opt/OpenBLAS/include -I/project-files/include -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran -L/project-files/usr/lib64 -lRmath -L/usr/lib -lR


nngp: nngp.o kvpar.o
	$(CXX) $(OPTFLAGS) -o nngp nngp.o kvpar.o $(LIBS) 

nngp.o: nngp.cpp
	$(CXX) $(OPTFLAGS) nngp.cpp -c $(LIBS)

kvpar.o: ../libs/kvpar.cpp
	$(CXX) $(OPTFLAGS) ../libs/kvpar.cpp -c

clean:
	rm -rf *o nngp
