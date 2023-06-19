CC = g++
LIBS = 

INC = -I/media/cubist/files/all/cpp/tools/eigen/eigen-eigen-b3f3d4950030/Eigen \
			-I/media/cubist/files/all/cpp/tools/myheaders \

CFLAGS = $(INC)

DEPS = $(find /media/cubist/files/all/cpp/tools/myheaders -name '*.h')
COBJS = $(find /media/cubist/files/all/cpp/tools/myheaders -name '*.o') \
				qftLatticeFunctions.cc

%.o : %.cc %.cpp $(DEPS) qftLattice.o
	$(CC) -c `root-config --libs --cflags` -o $@ $< $(CFLAGS)

qftLattice: $(COBJS) qftLattice.cc
	$(CC) $^ -o $@ `root-config --libs --cflags` $(CFLAGS)

clean:
	rm *~ *.o
