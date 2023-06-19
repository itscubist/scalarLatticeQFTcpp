/********
*Author: Baran Bodur
*Date: 2019-03-23
*Description: Lattice Monte Carlo for Scalar Boson QFT with quartic interaction in 2-3-4 space-time
* dimensions. Ising model like single cluster algorithm for the sign of the field at each location.
* An additional update for the magnitude of the field at each lattice site. Periodic boundary
* conditions is applied. 
* 
* There are 3 correlation functions calculated:
*		1) sigmaVar - Variance of every point
*		2) chiPairCor - Pair correlation
*		3) fPairCorTime -cos(2*\pi*delta(time)/Length) weighted pair correlation
*
* Inputs are:
*  1) LATLEN: Cubic Lattice Size in Every Direction
*  2) DIM : Dimension (ONLY 2 or 3 or 4) of the lattice
*  3) KAPPA : (2*DIM+ma^2)^-1  or (depends on the first option) ALPHA: mass*latSpace^2 (ma^2)
*  4) COUPLING: coupling constant of the quartic interaction
*  5) NSWEEP: Number of sweeps
*  6) FILEVER: Version No of the Output File
*  7) Options (will be a string of 1s and 0s, Ex: 101) - Default is 0  
*  		- Option 1: Chooses KAPPA or ALPHA as the third input parameter
*  			0) KAPPA
*  			1) ALPHA
*  		- Option 2: Histogram/Graph Disable Option
*  			0) Graphing/Histograming enabled
*  			1) Graphing/Histograming disabled
*			- Option 3: Shut Up Mode (Only the final results will be printed when 1)
*				0) Certain Debug Steps Will Be Printed 
*				1) Only Final Output Will Be Printed
*			- Option 4: Write The Field Configuration to File (used to start from equilibrium)
*				0) Do not Write
*				1) Write to the same filename as .root file but with .field extension instead
*			- Option 5: Read Field Configuration From File 
*				0) Init. field randomly without any input file
*				1) Init. by reading the field configuration file (run and save dim. must be the same)
*		8) INFILE: Name of the File to Read if option 5 is set
*
*  The output will both print the immediate results and save certain histograms in a .root file
*  for inspection/debug purposes 
*  (Filename: L(LATLEN)_D(DIM)_K(KAPPA)_G(COUPLING)_N(NSWEEP)_V(FILEVER).root or .field)
********/
// c++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <queue>
#include <bitset>
// ROOT headers
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"

// Eigen headers
#include "Eigen"

// Project Headers
#include "qftLattice.h"

using namespace std;
using namespace TMath;

//Main Function
int main(int argc, const char *argv[]) {
	//****** Obtain Inputs and Store in INPUT object
	if (argc<8) {
		cout <<"Usage: ./qftLattice LATLEN DIM KAPPA/ALPHA COUPLING NSWEEP FILEV OPTIONS INFILE"<< endl;
		return 0;
	} // end of usage printing statement
	QftInputPar INPUT; // define INPUT struct
	TString optString = (TString) argv[7]; // temp. option string
	for(unsigned int i=0;i<OPTLEN;i++) INPUT.OPTIONS[i] = (optString[i] == '1'); //init options array
	INPUT.LATLEN = (unsigned int) atoi(argv[1]); // lattice length
	INPUT.DIM = (unsigned int) atoi(argv[2]); INPUT.TWODIM = 2*INPUT.DIM; // dimension and 2*dim
	INPUT.LATSIZE = (unsigned int) pow(INPUT.LATLEN,INPUT.DIM);
	if(INPUT.OPTIONS[0] == true) { // init. kappa and alpha based on first option and relevant input
		INPUT.ALPHA = (double) atof(argv[3]);
		INPUT.KAPPA = 1.0/((double)(2*INPUT.DIM) + INPUT.ALPHA);
	}  
	else if(INPUT.OPTIONS[0] == false) {
		INPUT.KAPPA = (double) atof(argv[3]);
		INPUT.ALPHA = (double)1.0/INPUT.KAPPA - (double)(2*INPUT.DIM);
	} // end of init. kappa and alpha
	if(!INPUT.OPTIONS[2]) cout << "Kappa: " << INPUT.KAPPA << endl; 
	INPUT.COUPLING = (double) atof(argv[4]); // coupling const of the interaction term G
	INPUT.NSWEEP = (unsigned int) atoi(argv[5]); // number of sweep
	INPUT.FILEVER = (unsigned int) atoi(argv[6]); // saved file version
	INPUT.INFILE = (TString) argv[8]; // input file name
	INPUT.OUTFILE[0] = Form("outRootFiles/L%d_D%d_K%.2f_G%.2f_N%d_V%d.root", // root output filename
			INPUT.LATLEN,INPUT.DIM,INPUT.KAPPA,INPUT.COUPLING,INPUT.NSWEEP,INPUT.FILEVER);
	INPUT.OUTFILE[1] = Form("outFieldFiles/L%d_D%d_K%.2f_G%.2f_N%d_V%d.field", // field output filename
			INPUT.LATLEN,INPUT.DIM,INPUT.KAPPA,INPUT.COUPLING,INPUT.NSWEEP,INPUT.FILEVER);

	if(!INPUT.OPTIONS[2]) cout << "Inputs are recorded.\nInitializing Lattice." << endl;	
	
	//******* Initializing lattice based on inputs
	TRandom3 RNG; RNG.SetSeed(0); // init global RNG with random seed
	Lattice fLattice(INPUT); // Creat the lattice
	if(!INPUT.OPTIONS[2]) cout << "Lattice Initialized." << endl;		
	/* for (int i = 0; i < INPUT.LATSIZE; i++) {
		cout << "Abs Field Value: " << fLattice.latArr[i]->fVal << endl;
		cout << "Field cos(2*pi*Time/LATLEN): " << fLattice.latArr[i]->trigTime[0] << endl;
		cout << "Field sin(2*pi*Time/LATLEN): " << fLattice.latArr[i]->trigTime[1] << endl;
		for (int nCtr = 0; nCtr < INPUT.TWODIM; nCtr++) {
			cout << "Abs Field Value of Neighbour: " << nCtr <<" "<< 
				fLattice.latArr[i]->nbrs[nCtr]->fVal << endl;
		}
	} */
	
	double stats[2][3]; // mean and std of observables
	double aveClSize=0;
	unsigned int updateCount = INPUT.LATSIZE/2; // updates per sweep
	unsigned int updateRatio = 5*(1+INPUT.DIM-2); //sign to value update ratio 
	// pass that many sweeps before affecting the stats
	unsigned int passCount = 1; 
	if(INPUT.LATSIZE/updateCount>1 && INPUT.OPTIONS[4]==0) passCount=INPUT.LATSIZE/updateCount;
	if(!INPUT.OPTIONS[2]) {
		cout << "Pass Count: " << passCount << endl;
		cout << "Update Count: " << updateCount << endl;
	}
	double observables[INPUT.NSWEEP+passCount][3]; // observables
	// Graphs and Saving	
	TFile* outFile = new TFile(INPUT.OUTFILE[0],"RECREATE"); outFile->cd(); // Create out root file
	TGraph* obsGraph[4];
	for (int i = 0; i < 4; i++) obsGraph[i] = new TGraph(INPUT.NSWEEP+passCount);
	
	
	for (unsigned int swpCtr = 0; swpCtr < INPUT.NSWEEP+passCount; swpCtr++) {
		double tempObs[3] = {0,0,0};
		double clSize=0;
		for(unsigned int upCtr = 0; upCtr < updateCount; upCtr++) {
			fLattice.signUpdate(INPUT, RNG);
			//cout << "Done Sign Update: " << upCtr+swpCtr*updateCount << endl;
			if(upCtr % updateRatio == 0) fLattice.valUpdate(INPUT, RNG);
			//fLattice.ineffMeas(INPUT, RNG);
			//cout << "Done Value Update: "<< upCtr+swpCtr*updateCount <<  endl;
			clSize+=(double)fLattice.cSize;
			for (unsigned int i=0; i <3;i++) tempObs[i]+=fLattice.observable[i];
		} // end of in sweep loop
		for (unsigned int i=0; i <3;i++) 
			observables[swpCtr][i]=tempObs[i]/((double)(updateCount));
		clSize/=(double)updateCount;
		if(!INPUT.OPTIONS[2]) { // printing condition
			//if(swpCtr%100==0) 
			cout << swpCtr << " sweeps performed!" << endl;		
			cout << "Cluster Size: " << clSize << endl;		
			cout << "Sigma: " << observables[swpCtr][0] << endl;		
			cout << "Chi: " << observables[swpCtr][1] << endl;		
			cout << "F: " << observables[swpCtr][2] << endl;		
		} // end of printing condtion
		
		if(swpCtr>=passCount) { // If enough steps passed use info in stats calculation
			aveClSize+=(double)clSize;
			for(int i=0;i<3;i++) { // For stats
				stats[0][i]+=observables[swpCtr][i];
				stats[1][i]+=pow(observables[swpCtr][i],2);
			} // end of stats
		} // end of stats if
		// Nextt 2 Lines For graphing
		for(int i=0;i<3;i++) obsGraph[i]->SetPoint(swpCtr,(double)swpCtr,observables[swpCtr][i]);
		obsGraph[3]->SetPoint(swpCtr,(double)swpCtr,clSize);
	} // end of all sweeps loop

	aveClSize/=((double)(INPUT.NSWEEP));
	for(int i=0;i<3;i++) { // For stats
		stats[0][i]/=((double)(INPUT.NSWEEP));
		stats[1][i]/=((double)(INPUT.NSWEEP));
		stats[1][i]-=pow(stats[0][i],2);
		stats[1][i]=sqrt(stats[1][i]/((double)(INPUT.NSWEEP)));	
	} // end of stats calculations

	// Calculate M(L) and its error
	double massL[2]; 
	double tempChiF = 1/sqrt(stats[0][1]/stats[0][2]-1.0);
	massL[0] = 2*sin(Pi()/((double)INPUT.LATLEN))*tempChiF;
	massL[1] = sin(Pi()/((double)INPUT.LATLEN))*pow(tempChiF,3)*sqrt( pow(stats[1][1]/stats[0][2],2)
			+ pow((stats[1][2]*stats[0][1]/pow(stats[0][2],2)),2) );

	if(!INPUT.OPTIONS[2]) { // printing condition
		cout << "Ave. Cluster. Size: " << aveClSize << endl;		
		cout << "Sigma: " << stats[0][0] << " +- " << stats[1][0] << endl;		
		cout << "Chi: " << stats[0][1] << " +- " << stats[1][1] << endl;		
		cout << "F: " << stats[0][2] << " +- " << stats[1][2] << endl;
		cout << "M(L): " << massL[0] << " +- " << massL[1] << endl;
	} // end of printing condtion

	for(int i=0;i<4;i++) {
		obsGraph[i]->Draw("A*");
		obsGraph[i]->SetLineColor(4);
		obsGraph[i]->Write();
	} // end of graph writing
	outFile->Close(); // close root file	
	if(INPUT.OPTIONS[3]) fLattice.writeField(INPUT); // Save field for possible later use	
	return 0;
} // end of main
