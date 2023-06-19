/********
*Author: Baran Bodur
*Date: 2018-11-06
*Description: Monte Carlo of 3D Ising model a with single cluster algorithm. The observable we are
* looking at is spin susceptibility. Inputs are: J (the energy term per spin-spin in
* hamiltonian) multiplied by beta of interest (jBeta); length of the lattice in 1D (assume
* cubic lattice with periodic boundary conditions); number of sweeps (how many times an
* observable will be collected) - Notice in single cluster algorithm average of certain number
* of measurements on single clusters are called an observable, so total number of spin flips
* will be more than the sweep entered...
*
*  The output will both print the immediate results and give certain histograms in a .root file
********/
// c++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include<queue>
#include <bitset>
// ROOT headers
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <TRandom2.h>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"

using namespace std;
using namespace TMath;

const unsigned int MAXEXACT = 2; // Max Length for Exact Calculation
const unsigned int MAXBITS = pow(MAXEXACT,3); //Max bits in bitset, can do only 2 by 2 by 2 for now...
const unsigned int NNBRS = 6; //number of neighbours
struct SpinStruct;
//Spin Struct for Easily Accessing Values, Flags, Neighbours
struct SpinStruct{
	bool value=0; //Spin Value 0 or 1
	bool flag=0; //Spin in Cluster or No
	SpinStruct *nbrs[NNBRS]={NULL,NULL,NULL,NULL,NULL,NULL}; //Pointer to Neighbours
	unsigned int xVal,yVal,zVal;
};
//Function Declarations
inline void addN2Queue(queue<SpinStruct*> &Q, SpinStruct &curSpin, double pBond, TRandom3* inRand, 
		bool spinVal, bool verbose);
void initSpins(vector<vector<vector<SpinStruct>>> &spinsIn, const unsigned int length, bool isRand); 
long double exactSus(double jBeta, unsigned int LATLEN, TH1D &hist1, TGraph &graph1, double g1norm);
long double setSpins(vector<vector<vector<SpinStruct>>> &spins,
		unsigned int LATLEN,unsigned int stateID);	
void printSpins(vector<vector<vector<SpinStruct>>> &spins,unsigned int LATLEN);
unsigned int getStateID(vector<vector<vector<SpinStruct>>> &spins,unsigned int LATLEN);

//Main Function
int main(int argc, const char *argv[]) {
	// Obtain Input
	if (argc<4) {
			cout << "Usage: ./ising3D jBeta LATLEN NSWEEP (optionally save file version)" << endl;
	}
	double jBeta = (double)atof(argv[1]);
	const unsigned int LATLEN = (unsigned short) atoi(argv[2]);
	const unsigned int NSWEEP = (unsigned int) atoi(argv[3]);
	unsigned short version = 0;
	if(argc >= 5) version = atoi(argv[4]);
	unsigned int susDistBins = pow(LATLEN,3)*100;
	if(argc>=6) susDistBins = atoi(argv[5]);
	bool fitOut = false;
	if(argc>=7) fitOut = atoi(argv[6]);
	//Calculate these quantities only once
	double pBond = 1.000 - Exp(-2.0*jBeta); // probability of putting a bond
	if(!fitOut) cout << "Probability of Bond is: " << pBond << endl;
	unsigned int Volume = pow(LATLEN,3); // Volume
	//Conditions for measurement first one is conditioning total number of updates needed for a
	//measurement, second one is how many of them we skip after before starting a new measurement
	// 3rd input is the number of updates in a sweep
  unsigned int swpCnd[3] = {0,0,100}; //Based on distributions
  unsigned int printLim = 0; // How many sweeps are printed out
	if(fitOut) printLim = 0; // No prints if this is a fitting run
	//ROOT Plotting Inits
	TString fileName = Form("./rootfiles/ising3D_jb%f_L%d_S%d_v%d.root",jBeta,LATLEN,NSWEEP,version);
	TFile outFile(fileName,"RECREATE");
	//TH1::SetDefaultSumw2();
	TH1D susDist("susDist","Susceptibility Distribution",susDistBins,0,Volume+1);
	TH1D SijNdist("SijNdist","Density of States for 2x2x2 Ising Model",50,-25,25);
	TH1I withinSweep("withinSweep","Distribution of Cluster Sizes",Volume+1,0,Volume+1);
	TH1I withinSweepSel("withinSweepSel","Distribution of Cluster Sizes",Volume+1,0,Volume+1);
	TH1I statesMC("statesMC","MC Distribution of States",pow(2,Volume)+1,0,pow(2,Volume)+1);
	TH1I checkDist("checkDist","Dist of Queue Size",7,0,7);
	TH1I randCheck("randCheck","Random spin selection check",8,0,8);
	TGraph statesExact(pow(2,Volume));
  //Initialize spin lattice randomly, the last index:if 0 spin,if 1 in cluster flag
	vector<vector<vector<SpinStruct>>> spins
		(LATLEN,vector<vector<SpinStruct>>(LATLEN,vector<SpinStruct>(LATLEN)));
	initSpins(spins,LATLEN,true);
	if(!fitOut) cout << "Initialized Spins" << endl;
	//Create random number generator
	TRandom3 randSpin;
	randSpin.SetSeed(0);
	//Define within update variables
	SpinStruct *curSpin;
	//Initialize cluster and temp queues
	double clSize; // cluster size within an update
	double sumClSize; // Sum of all cluster sizes obtained within the sweep
	double clSizeMsCnt; // Number of cluster size obtained in sweep
	queue<SpinStruct*> clusterQ; // Keeps the spins in cluster
	queue<SpinStruct*> tempQ; // Keeps neighbours to check
	//Initialize Observable Array
	bool measFlag;
	double NMEAS = (double)NSWEEP;
	double spinSus[(int)NMEAS];
	double spinSusMean=0,spinSusStd=0,spinSusRMS=0,spinSusError;
	bool spinVal;
	//Main Loop (Sweeps)
	for(unsigned int swpctr=0;swpctr<NSWEEP;swpctr++) {
		//cout << "Starting Sweep: " << swpctr+1 << endl;
		sumClSize = 0;
		clSizeMsCnt = 0;
		measFlag = true;
		while(true) { // While loop for updates within a sweep
			// ************ FOR DEBUG ONLY, do not forget to remove!!!!
			//initSpins(spins,LATLEN,false);
			//cout << "Starting Update:" << clSizeMsCnt+1 << " in Sweep: " << swpctr+1 << endl;
			curSpin=&spins[randSpin.Integer(LATLEN)][randSpin.Integer(LATLEN)][randSpin.Integer(LATLEN)];
			randCheck.Fill(4*curSpin->xVal+2*curSpin->yVal+curSpin->zVal);
			//Add neighbouring spins to queue if same spin and if bonded	
			spinVal = curSpin->value;
			curSpin->flag = true; //Set flag of first spin
			if(swpctr< printLim && LATLEN<=MAXEXACT) printSpins(spins,LATLEN);
			addN2Queue(tempQ,*curSpin,pBond,&randSpin,spinVal,(swpctr<1?printLim:0));
			if(swpctr < printLim && LATLEN<=MAXEXACT) cout << "Queue Size: " << tempQ.size() << endl;
			checkDist.Fill(tempQ.size());
			clusterQ.push(curSpin);
			clSize = 1.0;
			//cout << "Queue Size: " << tempQ.size() << endl;
			//curSpin->value = !curSpin->value; //Flip the spin
			while(!tempQ.empty()) { //Keep expanding the cluster until the queue is empty (nowhere to look)
				curSpin = tempQ.front(); // Get the first element of the queue
				tempQ.pop(); //Remove the element obtained from queue
				if(swpctr < printLim && LATLEN<=MAXEXACT) printSpins(spins,LATLEN);
				addN2Queue(tempQ,*curSpin,pBond,&randSpin,spinVal,(swpctr<1?printLim:0)); //
				if(swpctr < printLim && LATLEN<=MAXEXACT) cout << "Queue Size: " << tempQ.size() << endl;
				clusterQ.push(curSpin);
				clSize+=1.0;
			}
			sumClSize+=clSize; // Increase summation
			clSizeMsCnt+=1.0;	
			
			if(swpctr < printLim && LATLEN <=MAXEXACT) {
				cout << "Cl Size of the above chain of events is: " << clSize << endl;
			}
			
			if(swpctr < printLim && LATLEN<=MAXEXACT) {
				printSpins(spins,LATLEN);
				cout << "And after flipping and unmasking we have the following config: " << endl;
			}
			while (!clusterQ.empty()) {
				curSpin = clusterQ.front();
				clusterQ.pop();
				curSpin->value = !spinVal;
				curSpin->flag = false;
			}
			
			if(swpctr < printLim && LATLEN<=MAXEXACT) printSpins(spins,LATLEN);
			withinSweep.Fill(clSize);	
			
			// The loop that determines how many updates to skip
			if(!measFlag && sumClSize > swpCnd[1]) {
				measFlag = true;
				sumClSize = 0;
				clSizeMsCnt = 0;
			} else if(measFlag) {
					withinSweepSel.Fill(clSize);
					if(clSizeMsCnt > swpCnd[2]) break;
					//if(sumClSize > swpCnd[0]) break;	
			}
		}
		if(LATLEN<=MAXEXACT) {
			statesMC.Fill(static_cast<int>(getStateID(spins,LATLEN))); // Fill state
		}
		//Average of cluster size is spin susceptibility
		//cout << "Cluster Sum: " << sumClSize << " Meas Count: "  << clSizeMsCnt << endl;
		spinSus[swpctr] = sumClSize/clSizeMsCnt;
		susDist.Fill(spinSus[swpctr]);
		// cout << "Spin Sus: " << spinSus[swpctr] << endl; // For Debug
	}
	// Calculate Average/Std. Dev.
	for (unsigned int swpctr = 0; swpctr < NMEAS; swpctr++) {
		spinSusMean+=spinSus[swpctr];
		spinSusRMS+=pow(spinSus[swpctr],2);
	}
	spinSusMean/=NMEAS;
	spinSusRMS/=NMEAS;
	spinSusStd = sqrt(spinSusRMS - pow(spinSusMean,2));
	spinSusError = spinSusStd/sqrt(NMEAS);
  if(!fitOut) {
  	if(LATLEN <= MAXEXACT) { //Print Exact Result
  		long double spinSusExact = exactSus(jBeta,LATLEN,SijNdist,statesExact,statesMC.Integral()); 
  		cout << "Exact Result is: " << setprecision(5) << spinSusExact << endl;
		}
		cout << "spinSusMean: " << spinSusMean << endl;
		cout << "spinSusError: " << spinSusError << endl;
		cout << "1st State Still Has: " << statesMC.GetBinContent(1) << endl;
	//Distributions to compare for DEBUG
	//for (int i = 0; i <= NNBRS; i++) {
		//cout<<"p"<<i<<" Actual: "<<checkDist.Integral()*ROOT::Math::binomial_pdf(i,pBond,6)
		//		<<" Current: "<<checkDist.GetBinContent(i+1) << endl;
		//}
	}
	// If fitting run only print results in a row: 
	if(fitOut) cout << LATLEN << " " << spinSusMean << " " << spinSusError << endl;
	//Write To File
	outFile.cd();
	susDist.Write();
	SijNdist.Write();
	withinSweep.Write();
	withinSweepSel.Write();
	randCheck.Write();
	if(LATLEN <= MAXEXACT) {
		statesMC.GetXaxis()->SetTitle("State ID (0-256)");
		statesMC.SetLineWidth(2);
		statesMC.SetLineColor(4);
		statesMC.Write();
		statesExact.SetLineWidth(2);
		statesExact.SetLineColor(2);
		statesExact.Write();
	}
	checkDist.Write();
	outFile.Close();
	return 0;
}

//Function to decide add neighbour pointers to queue
void addN2Queue(queue<SpinStruct*> &Q, SpinStruct &curSpin, double pBond, TRandom3* inRand, 
		bool spinVal, bool verbose) {
	double tempRand;
	if(verbose){
		cout<<"Considering Neighbors of Spin: "<<curSpin.xVal<<" "<<curSpin.yVal<<" "<<curSpin.zVal<<endl;
	}
	for(unsigned int nbrcnt=0;nbrcnt<NNBRS;nbrcnt++) {
		if(verbose) {
			cout << "Checking Neighbour: " << nbrcnt+1 << endl;
			cout << "Location is: " << curSpin.nbrs[nbrcnt]->xVal << " " << curSpin.nbrs[nbrcnt]->yVal <<
					" " << curSpin.nbrs[nbrcnt]->zVal << endl;
		}
		if(curSpin.nbrs[nbrcnt]->value == spinVal) {
			if(verbose)
				cout << "Spins matched: " << curSpin.nbrs[nbrcnt]->value << " and: " << spinVal << endl;
			if(!curSpin.nbrs[nbrcnt]->flag)  {
				tempRand = inRand->Rndm();
				if(verbose) {
					cout << "Not added to cluster before" << endl;
					cout << "Rolled: " << tempRand << endl;
				}
				if(tempRand <= pBond) {
					if(verbose) cout << "Spin is added to the cluster and marked" << endl;
					curSpin.nbrs[nbrcnt]->flag=true;
					Q.push(curSpin.nbrs[nbrcnt]); //Store neighbour pointers in queue
				}
			}
			else if(verbose){
				cout << "Spin was already in the cluster! skipping..." << endl;
			}
		}
	}
}

// Function to initalize spins values randomly, and sets periodic neighbouring relations
void initSpins(vector<vector<vector<SpinStruct>>> &spinsIn,const unsigned int length, bool isRand) { 
	TRandom2 initializer;
	initializer.SetSeed(0);
	for (unsigned int x = 0; x < length; x++) {
		for (unsigned int y = 0; y < length; y++) {
			for (unsigned int z = 0; z < length; z++) {
				if(isRand) spinsIn[x][y][z].value = initializer.Integer(2);
				else spinsIn[x][y][z].value = 0;
				spinsIn[x][y][z].flag = 0;
				spinsIn[x][y][z].xVal = x;
				spinsIn[x][y][z].yVal = y;
				spinsIn[x][y][z].zVal = z;
				if(x == 0) {	 //Special Attention if x is 0
					spinsIn[x][y][z].nbrs[0] = &spinsIn[length-1][y][z];
					spinsIn[x][y][z].nbrs[1] = &spinsIn[x+1][y][z];
				}
				else if(x==(length-1)) { // Special Attention if x=length-1
					spinsIn[x][y][z].nbrs[0] = &spinsIn[x-1][y][z];
					spinsIn[x][y][z].nbrs[1] = &spinsIn[0][y][z];
				}
				else {//Default neighbours are spins at x+1 and x-1
					spinsIn[x][y][z].nbrs[0] = &spinsIn[x-1][y][z];
					spinsIn[x][y][z].nbrs[1] = &spinsIn[x+1][y][z];
				}
				if (y==0){
					spinsIn[x][y][z].nbrs[2] = &spinsIn[x][length-1][z];
					spinsIn[x][y][z].nbrs[3] = &spinsIn[x][y+1][z];
				}
				else if(y==(length-1)) {
					spinsIn[x][y][z].nbrs[2] = &spinsIn[x][y-1][z];
					spinsIn[x][y][z].nbrs[3] = &spinsIn[x][0][z];
				}
				else {
					spinsIn[x][y][z].nbrs[2] = &spinsIn[x][y-1][z];
					spinsIn[x][y][z].nbrs[3] = &spinsIn[x][y+1][z];
				}
				if(z==0) {
					spinsIn[x][y][z].nbrs[4] = &spinsIn[x][y][length-1];
					spinsIn[x][y][z].nbrs[5] = &spinsIn[x][y][z+1];
				}
				else if(z==(length-1)) {	
					spinsIn[x][y][z].nbrs[4] = &spinsIn[x][y][z-1];
					spinsIn[x][y][z].nbrs[5] = &spinsIn[x][y][0];
				}
				else {
					spinsIn[x][y][z].nbrs[4] = &spinsIn[x][y][z-1];
					spinsIn[x][y][z].nbrs[5] = &spinsIn[x][y][z+1];
				}
			}			
		}		
	}
}
//Exact calculation for 3D with input jBeta
long double exactSus(double jBeta, unsigned int LATLEN, TH1D &hist1, TGraph &graph1, double g1norm) {
	//Definitions
	ofstream exact("exacSol222.txt");
	unsigned int VOLUME = pow(LATLEN,3);
	unsigned int NSTATES = pow(2,VOLUME);
	vector<long double> probs(NSTATES);
	vector<long double> expecTerms(NSTATES);
	queue<SpinStruct*> tempQ;
	SpinStruct *curSpin;
	long double sumSijN; //sum of all neighbouring spins for a given configuration
	long double sumSij; //sum of all spin-spin multiplications
	long double partFunc=0,susExact=0,dummy,firstProb;
	// Declare & Initialize spin configuration
	vector<vector<vector<SpinStruct>>> spins
		(LATLEN,vector<vector<SpinStruct>>(LATLEN,vector<SpinStruct>(LATLEN)));
	initSpins(spins,LATLEN,false);
	
	cout << "Initialized Spins For Exact Solutions" << endl;
	
	// Loop over all possible configurations
	for (unsigned int stctr=0;stctr < NSTATES;stctr++) {
		sumSij = setSpins(spins,LATLEN,stctr); //set Spins to the new configuration 
		//cout << "Spins Set to Configuration: " << stctr+1 << endl;
		//cout << "Magnetization square per Volume of the state is: " << sumSij << endl;
		curSpin = &spins[0][0][0]; //Init current spin pointer
		sumSijN = 0; // Re-zero spin spin configurations

		for (unsigned int x = 0; x < LATLEN; x++) {
			for (unsigned int y = 0; y < LATLEN; y++) {
				for (unsigned int z = 0; z < LATLEN; z++) {
					for(unsigned int nbrcnt=0; nbrcnt<NNBRS; nbrcnt++) {
						if(spins[x][y][z].value == spins[x][y][z].nbrs[nbrcnt]->value) sumSijN+=0.5;
						else sumSijN-=0.5;
					}
				}
			}
		}
		//cout << "sumSijN of the state is: " << sumSijN << endl;
		// Output
		hist1.Fill(sumSijN);
		exact << stctr+1 << " " << sumSijN << " " << sumSij << endl;
		dummy = Exp(jBeta*sumSijN);
		probs[stctr] = dummy;
		expecTerms[stctr] = (long double)sumSij*dummy;
		
		//cout << "prob: " << setprecision(5) << probs[stctr] << endl;
		//cout << "expected: " << expecTerms[stctr] << endl;
	}
	// Sort the vectors to make calculation more precise
	//sort(probs.begin(),probs.begin()+NSTATES);
	//sort(expecTerms.begin(),expecTerms.begin()+NSTATES);	
	for (unsigned int stctr=0;stctr < NSTATES;stctr++) {
		partFunc+=probs[stctr];
		susExact+=expecTerms[stctr];
	}
	for (unsigned int stctr=0;stctr < NSTATES;stctr++) {
		graph1.SetPoint(static_cast<int>(stctr),(double)stctr,probs[stctr]*((double)g1norm/partFunc));
	}
	exact.close();
	return susExact/partFunc;
}

//Set spin configuration to binary equivalent of input number and set flags to 0
//Also returns the magnetization square per volume of the microstate
long double setSpins(vector<vector<vector<SpinStruct>>> &spins,unsigned int LATLEN,
		unsigned int stateID) {	
	long double VOLUME = pow((double)LATLEN,3);
	bitset<MAXBITS> spinBits(stateID);
	unsigned int bitCtr = 0;
	for (unsigned int x = 0; x < LATLEN; x++) {
		for (unsigned int y = 0; y < LATLEN; y++) {
			for (unsigned int z = 0; z < LATLEN; z++) {
				spins[x][y][z].value = spinBits[bitCtr];
				bitCtr++;
				spins[x][y][z].flag = false;
			}
		}
	}
	return pow((VOLUME-2*(long double)spinBits.count()),2)/VOLUME;
}

// Works for 2x2x2 for now...
void printSpins(vector<vector<vector<SpinStruct>>> &spins,unsigned int LATLEN) {
	int z = 0;
	cout << "******************" << endl;
	for (unsigned int x = 0; x < LATLEN; x++) {
		for (unsigned int y = 0; y < LATLEN; y++) {
			//for (unsigned int z = 0; z < LATLEN; z++) {
			cout << spins[x][y][z].value << " | " << spins[x][y][z+1].value;
			cout << "   |||||   " << spins[x][y][z].flag << " | " << spins[x][y][z+1].flag;
			//}
			cout << endl;
		}
	cout << "------------------" << endl;
	}
	cout << "******************" << endl;
}

unsigned int getStateID(vector<vector<vector<SpinStruct>>> &spins,unsigned int LATLEN) {
	bitset<MAXBITS> spinBits(0);
	unsigned int bitCtr = 0;
	for (unsigned int x = 0; x < LATLEN; x++) {
		for (unsigned int y = 0; y < LATLEN; y++) {
			for (unsigned int z = 0; z < LATLEN; z++) {
				spinBits.set(bitCtr,spins[x][y][z].value);
				bitCtr++;
			}
		}
	}
	return spinBits.to_ulong();
}
