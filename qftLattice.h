#ifndef QFTLATTICE_H	
	#define QFTLATTICE_H

// c++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
// ROOT headers
#include "TMath.h"
#include "TString.h"
#include "TRandom.h"
#include "TRandom3.h"
// Eigen headers
#include "Eigen"

const unsigned int OPTLEN=5; //length of options input


// Declare structs/classes
struct QftInputPar;
struct LatticePts;
class Lattice;

struct QftInputPar{ // INPUT Struct Definition
	unsigned int LATLEN;
	unsigned int DIM;
	unsigned int TWODIM;
	unsigned int LATSIZE;
	double KAPPA;
	double ALPHA;
	double COUPLING;
	unsigned int FILEVER;
	unsigned int NSWEEP;
	TString INFILE;
	TString OUTFILE[2];
	bool OPTIONS[OPTLEN];
};

struct LatticePts{ // Struct for a point in lattice, time is always 1st dimension
	bool fSign; 
	double fVal;
	double trigTime[2]; // 0: cosine, 1: sine of  (2*\pi*\timeDim)/LATLEN, cos(a-b)=cosa*cosb-sina*sinb
	LatticePts* nbrs[8] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL}; //pointers to neighbours
};

class Lattice{
public:
	// Public Functions
	Lattice(QftInputPar INPUT); // Constructor
	void writeField(QftInputPar INPUT); // Write field to input text file
	void signUpdate(QftInputPar INPUT, TRandom3& inRand); // Perform Sign Update of the Field
	void valUpdate(QftInputPar INPUT,TRandom3& inRand); // Perform Value Uppdate of the Field
	void ineffMeas(QftInputPar INPUT,TRandom3& inRand); // Inefficient Measure to Check
	// Variables
	double observable[3]; // Filled after running measObs with the current configs observables
	int cSize; // size of the clusters
	std::vector<LatticePts*> latArr; // Holds pointers to lattice points in 1D structure
private:
	// Private Functions	
	void initNeighbours(QftInputPar INPUT); // Neighbour init. and fill latArr
	// Actual 2,3,4 D vectors to init. neighbouring, should not be used after initialization
	std::vector<std::vector<LatticePts> > tempLattice2D;
	std::vector<std::vector<std::vector<LatticePts> > > tempLattice3D;
	std::vector<std::vector<std::vector<std::vector<LatticePts> > > > tempLattice4D;
};



#endif
