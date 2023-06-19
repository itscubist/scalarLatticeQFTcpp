using namespace std;

void qftFitter() {
	gStyle->SetOptFit(1111);
	// fit parameters and function
	const int nFit = 7;	
	double alpha[nFit] = {-0.55,-0.56,-0.57,-0.58,-0.59,-0.60,-0.61};
	double massL[nFit] = {0.202,0.186,0.170,0.153,0.136,0.118,0.100};
	double massLErr[nFit] = {0.001,0.001,0.001,0.001,0.001,0.001,0.001};
	//double alpha[nFit] = {-0.55,-0.60,-0.65,-0.70,-0.75};
	//double massL[nFit] = {0.578,0.531,0.478,0.421,0.354};
	//double massLErr[nFit] = {0.002,0.002,0.002,0.001,0.001};
	TF1 *fitFunc = new TF1("fitFunc","[0]*(x-[1])**[2]",-0.50,-0.80);
	fitFunc->SetParNames("f_{0}","#alpha_{c}","#nu");
  fitFunc->SetParameter(0,2);
	fitFunc->SetParameter(1,-0.85);
	fitFunc->SetParameter(2,1);
	//fitFunc->FixParameter(2,1);
	
	//Look at the fit function
	TCanvas* canny0 = new TCanvas("canny0","canny0",800,600);
	canny0->cd();
	fitFunc->Draw();

	// create fit graph
	TGraphErrors* fitGraph = new TGraphErrors(nFit,alpha,massL,NULL,massLErr);
	fitGraph->Fit("fitFunc");
	
	TCanvas* canny1 = new TCanvas("canny1","canny1",800,600);
	canny1->cd();
	fitGraph->SetTitle("d=3, L=12, g=0.01, M(L) vs. #alpha");
	fitGraph->GetXaxis()->SetTitle("#alpha");
	fitGraph->GetYaxis()->SetTitle("M(L)");
	fitGraph->SetMarkerSize(1);
	fitGraph->SetMarkerStyle(20);
	fitGraph->Draw("AP");
	canny1->Draw();

	// For spontaneous symmetry breaking part
	const int nL = 4;
	double lengths[nL] = {16,24,32,48};
	double m0p7Chis[nL] = {195,219,225,226};
	double m0p7ChisErr[nL] = {3,2,1,1};
	double m0p8Chis[nL] = {3144,10608,25029,84335};
	double m0p8ChisErr[nL] = {11,15,43,67};
	TGraphErrors* m0p7Graph = new TGraphErrors(nL,lengths,m0p7Chis,NULL,m0p7ChisErr);
	TGraphErrors* m0p8Graph = new TGraphErrors(nL,lengths,m0p8Chis,NULL,m0p8ChisErr);
	TMultiGraph* multi = new TMultiGraph();

	
	TF1 *fitFuncL = new TF1("fitFuncL","[0]*x**[1]",15,50);
	fitFuncL->SetParNames("coeff","power");
  fitFuncL->SetParameter(0,100);
	fitFuncL->SetParameter(1,3);
	m0p8Graph->Fit("fitFuncL");
	m0p7Graph->Fit("fitFuncL");

	m0p7Graph->SetMarkerStyle(20);
	m0p7Graph->SetMarkerColor(4); // blue
	m0p7Graph->SetMarkerSize(1.5);
	m0p7Graph->SetLineWidth(2);
	m0p7Graph->SetLineColor(4);
	m0p8Graph->SetMarkerStyle(21);
	m0p8Graph->SetMarkerColor(2); // red
	m0p8Graph->SetMarkerSize(1.5);
	m0p8Graph->SetLineWidth(2);
	m0p8Graph->SetLineColor(2);

	multi->SetTitle("d=3, g=0.01, L dependence at different #alpha; Lat. Len (L); #chi");

	multi->Add(m0p7Graph);
	multi->Add(m0p8Graph);

	TCanvas* canny2 = new TCanvas("canny2","canny2",800,600);
	canny2->SetLogy();
	canny2->cd();
	multi->Draw("AP");
	canny2->Draw();
}
