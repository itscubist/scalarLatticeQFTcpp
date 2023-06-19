/********
*Author: Baran Bodur
*Date: 2019-03-24
*Description: Definition of Necessary Functions For qftLattice
*
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
#include "TRandom.h"
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


// Lattice Class Constructor
Lattice::Lattice(QftInputPar INPUT) {
	latArr.reserve(INPUT.LATSIZE); // Reserve Memory for All Lattice Points
	initNeighbours(INPUT); // Arrange neighbours and time info
	if(INPUT.OPTIONS[2]==false) cout<<"Neighbours Initialized"<<endl;
	if(INPUT.OPTIONS[4]==false) { // Set initial fields randomly
		TRandom3 initRand; initRand.SetSeed(0); // Create a rng with random seed
		double tempRand; // temp variable
		for (int lPtr = 0; lPtr < INPUT.LATSIZE; lPtr++) { //
			tempRand = initRand.Gaus(0,10); // gaussion random var. with 0 mean sigma=100
			latArr[lPtr]->fSign = (tempRand > 0); // sign of field at that lattice ptr.
			latArr[lPtr]->fVal = Abs(tempRand); // absolute val of field at that lattice ptr.
		 } // end of loop over points for random init.
	}	// end of random init. condition
	else if(INPUT.OPTIONS[4] == true) { // init. from file
		ifstream inField(INPUT.INFILE); //create input stream
		for (int lPtr = 0; lPtr < INPUT.LATSIZE; lPtr++) { // loop over lattice points
			inField >> latArr[lPtr]->fSign >> latArr[lPtr]->fVal;	// read lattice pts
		} // end of loop over latttice points
		inField.close(); // close stream
	} // end of init. from file
	if(INPUT.OPTIONS[2]==false) cout<<"Field Lattice Initialized"<<endl;
} //end of contructor

// Lattice Class Neighbour (and Time) Initialiation Functions for 2,3 and 4D respectively

void Lattice::initNeighbours(QftInputPar INPUT) {
	double tempTrigTime[2];
	if(INPUT.DIM == 2) {
		tempLattice2D.assign(INPUT.LATLEN,std::vector<LatticePts>(INPUT.LATLEN) ); // init 2D vector
		for (unsigned int x = 0; x < INPUT.LATLEN; x++) {
			tempTrigTime[0] = cos(2*Pi()*((double)x/(double)INPUT.LATLEN));
			tempTrigTime[1] = sin(2*Pi()*((double)x/(double)INPUT.LATLEN));
			for (unsigned int y = 0; y < INPUT.LATLEN; y++) {
				tempLattice2D[x][y].trigTime[0] = tempTrigTime[0];
				tempLattice2D[x][y].trigTime[1] = tempTrigTime[1];	
				if(x == 0) {	 //Special Attention if x is 0
					tempLattice2D[x][y].nbrs[0] = &tempLattice2D[INPUT.LATLEN-1][y];
					tempLattice2D[x][y].nbrs[1] = &tempLattice2D[x+1][y];
				}
				else if(x==(INPUT.LATLEN-1)) { // Special Attention if x=INPUT.LATLEN-1
					tempLattice2D[x][y].nbrs[0] = &tempLattice2D[x-1][y];
					tempLattice2D[x][y].nbrs[1] = &tempLattice2D[0][y];
				}
				else {//Default neighbours are spins at x+1 and x-1
					tempLattice2D[x][y].nbrs[0] = &tempLattice2D[x-1][y];
					tempLattice2D[x][y].nbrs[1] = &tempLattice2D[x+1][y];
				} // end of 1st dim (time dim)
				if (y==0){ //2nd dim
					tempLattice2D[x][y].nbrs[2] = &tempLattice2D[x][INPUT.LATLEN-1];
					tempLattice2D[x][y].nbrs[3] = &tempLattice2D[x][y+1];
				}
				else if(y==(INPUT.LATLEN-1)) {
					tempLattice2D[x][y].nbrs[2] = &tempLattice2D[x][y-1];
					tempLattice2D[x][y].nbrs[3] = &tempLattice2D[x][0];
				}
				else {
					tempLattice2D[x][y].nbrs[2] = &tempLattice2D[x][y-1];
					tempLattice2D[x][y].nbrs[3] = &tempLattice2D[x][y+1];
				} //end of 2nd dim
				latArr.push_back(&tempLattice2D[x][y]);
			} // end x (time)
		} // end y
	} // end of 2D
	else if(INPUT.DIM == 3) {
		tempLattice3D.assign(INPUT.LATLEN,std::vector<std::vector<LatticePts> > (INPUT.LATLEN,
				std::vector<LatticePts>(INPUT.LATLEN) ) ); // 3D init
		for (unsigned int x = 0; x < INPUT.LATLEN; x++) {
			tempTrigTime[0] = cos(2*Pi()*((double)x/(double)INPUT.LATLEN));
			tempTrigTime[1] = sin(2*Pi()*((double)x/(double)INPUT.LATLEN));
			for (unsigned int y = 0; y < INPUT.LATLEN; y++) {
				for (unsigned int z = 0; z < INPUT.LATLEN; z++) {
					tempLattice3D[x][y][z].trigTime[0] = tempTrigTime[0];
					tempLattice3D[x][y][z].trigTime[1] = tempTrigTime[1];	
					if(x == 0) {	 //Special Attention if x is 0
						tempLattice3D[x][y][z].nbrs[0] = &tempLattice3D[INPUT.LATLEN-1][y][z];
						tempLattice3D[x][y][z].nbrs[1] = &tempLattice3D[x+1][y][z];
					}
					else if(x==(INPUT.LATLEN-1)) { // Special Attention if x=INPUT.LATLEN-1
						tempLattice3D[x][y][z].nbrs[0] = &tempLattice3D[x-1][y][z];
						tempLattice3D[x][y][z].nbrs[1] = &tempLattice3D[0][y][z];
					}
					else {//Default neighbours are spins at x+1 and x-1
						tempLattice3D[x][y][z].nbrs[0] = &tempLattice3D[x-1][y][z];
						tempLattice3D[x][y][z].nbrs[1] = &tempLattice3D[x+1][y][z];
					} // end of 1st dim (time dim)
					if (y==0){ //2nd dim
						tempLattice3D[x][y][z].nbrs[2] = &tempLattice3D[x][INPUT.LATLEN-1][z];
						tempLattice3D[x][y][z].nbrs[3] = &tempLattice3D[x][y+1][z];
					}
					else if(y==(INPUT.LATLEN-1)) {
						tempLattice3D[x][y][z].nbrs[2] = &tempLattice3D[x][y-1][z];
						tempLattice3D[x][y][z].nbrs[3] = &tempLattice3D[x][0][z];
					}
					else {
						tempLattice3D[x][y][z].nbrs[2] = &tempLattice3D[x][y-1][z];
						tempLattice3D[x][y][z].nbrs[3] = &tempLattice3D[x][y+1][z];
					} //end of 2nd dim
					if(z==0) { //3rd dim
						tempLattice3D[x][y][z].nbrs[4] = &tempLattice3D[x][y][INPUT.LATLEN-1];
						tempLattice3D[x][y][z].nbrs[5] = &tempLattice3D[x][y][z+1];
					}
					else if(z==(INPUT.LATLEN-1)) {	
						tempLattice3D[x][y][z].nbrs[4] = &tempLattice3D[x][y][z-1];
						tempLattice3D[x][y][z].nbrs[5] = &tempLattice3D[x][y][0];
					}
					else {
						tempLattice3D[x][y][z].nbrs[4] = &tempLattice3D[x][y][z-1];
						tempLattice3D[x][y][z].nbrs[5] = &tempLattice3D[x][y][z+1];
					} // end of 3rd dim
					latArr.push_back(&tempLattice3D[x][y][z]);
				} // end x (time)
			} // end y
		} // end z
	} //end of 3D
	else if(INPUT.DIM == 4) {
		tempLattice4D.assign(INPUT.LATLEN,std::vector<std::vector<std::vector<LatticePts> > > 
				(INPUT.LATLEN, std::vector<std::vector<LatticePts> > (INPUT.LATLEN,
				std::vector<LatticePts>(INPUT.LATLEN) ) ) ); // 4D init
		for (unsigned int x = 0; x < INPUT.LATLEN; x++) {
			tempTrigTime[0] = cos(2*Pi()*((double)x/(double)INPUT.LATLEN));
			tempTrigTime[1] = sin(2*Pi()*((double)x/(double)INPUT.LATLEN));
			for (unsigned int y = 0; y < INPUT.LATLEN; y++) {
				for (unsigned int z = 0; z < INPUT.LATLEN; z++) {
					for (unsigned int t = 0; t < INPUT.LATLEN; t++) {
						tempLattice4D[x][y][z][t].trigTime[0] = tempTrigTime[0];
						tempLattice4D[x][y][z][t].trigTime[1] = tempTrigTime[1];	
						if(x == 0) {	 //Special Attention if x is 0
							tempLattice4D[x][y][z][t].nbrs[0] = &tempLattice4D[INPUT.LATLEN-1][y][z][t];
							tempLattice4D[x][y][z][t].nbrs[1] = &tempLattice4D[x+1][y][z][t];
						}
						else if(x==(INPUT.LATLEN-1)) { // Special Attention if x=INPUT.LATLEN-1
							tempLattice4D[x][y][z][t].nbrs[0] = &tempLattice4D[x-1][y][z][t];
							tempLattice4D[x][y][z][t].nbrs[1] = &tempLattice4D[0][y][z][t];
						}
						else {//Default neighbours are spins at x+1 and x-1
							tempLattice4D[x][y][z][t].nbrs[0] = &tempLattice4D[x-1][y][z][t];
							tempLattice4D[x][y][z][t].nbrs[1] = &tempLattice4D[x+1][y][z][t];
						} // end of 1st dim (time dim)
						if (y==0){ //2nd dim
							tempLattice4D[x][y][z][t].nbrs[2] = &tempLattice4D[x][INPUT.LATLEN-1][z][t];
							tempLattice4D[x][y][z][t].nbrs[3] = &tempLattice4D[x][y+1][z][t];
						}
						else if(y==(INPUT.LATLEN-1)) {
							tempLattice4D[x][y][z][t].nbrs[2] = &tempLattice4D[x][y-1][z][t];
							tempLattice4D[x][y][z][t].nbrs[3] = &tempLattice4D[x][0][z][t];
						}
						else {
							tempLattice4D[x][y][z][t].nbrs[2] = &tempLattice4D[x][y-1][z][t];
							tempLattice4D[x][y][z][t].nbrs[3] = &tempLattice4D[x][y+1][z][t];
						} //end of 2nd dim
						if(z==0) { //3rd dim
							tempLattice4D[x][y][z][t].nbrs[4] = &tempLattice4D[x][y][INPUT.LATLEN-1][t];
							tempLattice4D[x][y][z][t].nbrs[5] = &tempLattice4D[x][y][z+1][t];
						}
						else if(z==(INPUT.LATLEN-1)) {	
							tempLattice4D[x][y][z][t].nbrs[4] = &tempLattice4D[x][y][z-1][t];
							tempLattice4D[x][y][z][t].nbrs[5] = &tempLattice4D[x][y][0][t];
						}
						else {
							tempLattice4D[x][y][z][t].nbrs[4] = &tempLattice4D[x][y][z-1][t];
							tempLattice4D[x][y][z][t].nbrs[5] = &tempLattice4D[x][y][z+1][t];
						} // end of 3rd dim
						if(t==0) { //3rd dim
							tempLattice4D[x][y][z][t].nbrs[6] = &tempLattice4D[x][y][z][INPUT.LATLEN-1];
							tempLattice4D[x][y][z][t].nbrs[7] = &tempLattice4D[x][y][z][t+1];
						}
						else if(t==(INPUT.LATLEN-1)) {	
							tempLattice4D[x][y][z][t].nbrs[6] = &tempLattice4D[x][y][z][t-1];
							tempLattice4D[x][y][z][t].nbrs[7] = &tempLattice4D[x][y][z][0];
						}
						else {
							tempLattice4D[x][y][z][t].nbrs[6] = &tempLattice4D[x][y][z][t-1];
							tempLattice4D[x][y][z][t].nbrs[7] = &tempLattice4D[x][y][z][t+1];
						} // end of 4th dim
					latArr.push_back(&tempLattice4D[x][y][z][t]);
					} // end x (time)
				} // end y
			} // end z
		} //end t	
	} // end of 4D
} // end of neighbour init. function

// Write field to input text file
void Lattice::writeField(QftInputPar INPUT) {
	ofstream outField(INPUT.OUTFILE[1]);
	for (int lCtr = 0; lCtr < INPUT.LATSIZE; lCtr++) {
		outField << latArr[lCtr]->fSign << " " << latArr[lCtr]->fVal << endl;
	}
	outField.close();
	if(!INPUT.OPTIONS[2]) cout << INPUT.LATSIZE << " field sign and values saved to file!" << endl;		
} // end of writeField
	
// Perform Sign Update of the Field
void Lattice::signUpdate(QftInputPar INPUT, TRandom3& inRand) { 
	LatticePts* curPts = latArr[inRand.Integer(INPUT.LATSIZE)]; // select a random site
	queue<LatticePts*> tempQ; // temporary queue to store lattice pts in queue
	double pNotBond; // probability of not adding a lattice point to queue
	bool tempSign = curPts->fSign; double tempVal = curPts->fVal; // save the current sign and value
	curPts->fSign = !tempSign; // flip the sign now to prevent visiting again
	tempQ.push(curPts); // push into the queue to add to cluster, and check neighbours
	
	double tempTrigTime[2] = {curPts->trigTime[0],curPts->trigTime[1]}; // save current time trigs
	observable[1] = 0; // contribution to observables from the 1st point in cluster
	observable[2] = 0; observable[0] = 0;
	cSize = 0;	
	while(!tempQ.empty()) { //Keep expanding the cluster until the queue is empty (nowhere to look)
		curPts = tempQ.front(); // Get the first element of the queue
		cSize++;
		observable[0] += pow(curPts->fVal,2);
		double tempMult = tempVal*curPts->fVal; // temporary storage, since will be used twice
		observable[1] += tempMult; // Calculate correlation observable
		observable[2] += tempMult*(tempTrigTime[0]*curPts->trigTime[0]+
				tempTrigTime[1]*curPts->trigTime[1]); // Calculate Time Weighted correlation observable
		tempQ.pop(); //Remove the element obtained from queue
		for(unsigned int nbrCtr=0; nbrCtr<INPUT.TWODIM; nbrCtr++) {
			//cout << "Checking Neighbour: " << nbrCtr+1 << endl;
			if(curPts->nbrs[nbrCtr]->fSign == tempSign) { // If same sign (otherwise not added to cluster)
				//cout << "Current fVal: " << curPts->fVal << endl;
				//cout << "Neighbour fVal: " << curPts->nbrs[nbrCtr]->fVal << endl;
				pNotBond = Exp(-2*INPUT.KAPPA*curPts->fVal*curPts->nbrs[nbrCtr]->fVal); // calculate prob
				//cout << "1-pBond: " << pNotBond << endl;
				if(inRand.Rndm() > pNotBond) { //If prob is small enough
					curPts->nbrs[nbrCtr]->fSign = !tempSign; // flip sign now, to prevent revisiting
					tempQ.push(curPts->nbrs[nbrCtr]); // add to queue, hence to the cluster
				} // end of lucky rng
			} // end of field sign equality check 
		} // end of neighbours of the current point loop
	} // end of queue emptying while loop

	//observable[1]/=(double)INPUT.LATSIZE; // calculate chi
	//observable[2]/=(double)INPUT.LATSIZE; // calculate f
	observable[0]/=(double)cSize;
} // end of signUpdate function

// Perform Value Uppdate of the Field
void Lattice::valUpdate(QftInputPar INPUT, TRandom3& inRand) {
	//observable[0] = 0;
	for(unsigned int lCtr = 0; lCtr<INPUT.LATSIZE; lCtr++) { // loop over lattice positions
		double fValNew = 0; //init neighbour contribution
		LatticePts* curPts = latArr[lCtr];
		for(unsigned int nbrCtr=0;nbrCtr<INPUT.TWODIM;nbrCtr++) {// loop over neighbours and add them
			fValNew+=curPts->nbrs[nbrCtr]->fVal*(curPts->nbrs[nbrCtr]->fSign==curPts->fSign ? 1.0:-1.0);
			//cout << "Neigh Val: " << curPts->nbrs[nbrCtr]->fVal << " Neigh Sign: " << 
			//		curPts->nbrs[nbrCtr]->fSign << endl;
		}
		fValNew*=INPUT.KAPPA; // multiply with KAPPA
		fValNew+=inRand.Gaus(0,1); // add normal rv
		//cout << "Proposed Value: " << fValNew << " Old Value: " << curPts->fVal << " Old Sign: " 
			//	<< curPts->fSign << endl;
		double fValSqr[2] = {pow(curPts->fVal,2),pow(fValNew,2)}; // old and new field squares
		double pAccept = Exp(INPUT.COUPLING*(pow(fValSqr[0],2)-pow(fValSqr[1],2))); 
		if(pAccept > inRand.Rndm()) {
			curPts->fVal = Abs(fValNew);
			curPts->fSign = (fValNew>0 ? curPts->fSign : !curPts->fSign);
			//observable[0] += fValSqr[1];
		} 
		//else observable[0] += fValSqr[0];
		//cout << "Resultant Value: " << curPts->fVal << " Res. Sign: " << curPts->fSign << endl;
	} // end of loop over lattice positions
	//observable[0]/=(double)INPUT.LATSIZE;
} // end of valUpdate function

// Inefficieny but Direct Measurement to Check
void Lattice::ineffMeas(QftInputPar INPUT, TRandom3& inRand) {
	observable[0] = 0; observable[1] = 0; observable[2]=0;
	for(unsigned int xCtr = 0; xCtr<INPUT.LATSIZE; xCtr++) { // loop over lattice positions
		LatticePts* xPts = latArr[xCtr];
		observable[0] += pow(xPts->fVal,2);
		for(unsigned int yCtr = 0; yCtr<INPUT.LATSIZE; yCtr++) {// loop over lattice positions
			LatticePts* yPts = latArr[yCtr];
			double tempMult = xPts->fVal*yPts->fVal*(xPts->fSign==yPts->fSign ? 1.0:-1.0);
			observable[1] += tempMult;
			observable[2] += tempMult*(xPts->trigTime[0]*yPts->trigTime[0]+
					xPts->trigTime[1]*yPts->trigTime[1]);
		}
	}
	for (int i = 0; i < 3; i++) observable[i]/=(double)INPUT.LATSIZE;	
}

