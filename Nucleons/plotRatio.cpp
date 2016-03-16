/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief       Plot ratio.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string>

#include "TROOT.h"
#include "TApplication.h"
#include "TComplex.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStopwatch.h"

///@{
char dataFile[] = "dataRatio.dat";           /**< Form factor data. */
///@}

char parametersFile[] = "parNucleons-FitD.dat";        ///< Input parameters.
char outputFile[] = "outNucleons-temp.dat";            ///< Output parameters.

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraph.cpp"
#include "modules/Nucleon3G.cpp"
#include "modules/NucleonFcn.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2016 Erik Bartoš" << std::endl;


	/// Start
	
	TApplication theApp("Ratios", &argc, argv);

	/// Plot
	std::cout << "\n> Plotting:" << std::endl;
	std::cout << "> Plotted parameters:          `" << parametersFile << "'" << std::endl;
	std::cout << "> Plotted form factor data:    `" << dataFile << "'" << std::endl;

	const int nPoints = 2500;
	double tMin;
	double tMax;
	double tStep;
	double tA;

	ExperimentalData W;
	W.ReadData(dataFile);
	int numG = W.size();
	std::cout << "\nNumber of plotted points:      " << numG << std::endl;
	std::vector<int> type = W.Type();
	std::vector<double> x = W.X();
	std::vector<double> val = W.Val();
	std::vector<double> errUp = W.ErrUp();
	std::vector<double> errDown = W.ErrDown();
	std::vector<double> theta = W.Theta();
	std::vector<double> energy = W.Energy();	

	FFactor pPlot(12);
	pPlot.LoadParameters(parametersFile);

	double plotRX[nPoints], plotRpY[nPoints];
	int numE = W.typeN(9);
	double dataX4[numE];
	double dataY4[numE];
	double dataU4[numE];
	double dataD4[numE];
	double tZ = 0.;
	
	tMin = 3.5;
	tMax = 10.0;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	TComplex t,hA,hB,hC,hD;

	int k4 = 0;
	for (int i = 0; i < W.size(); i++) {
		if (type[i] == 9) {
			dataX4[k4] = x[i];
			dataY4[k4] = val[i];
			dataD4[k4] = errDown[i];
			dataU4[k4] = errUp[i];
			++k4;
		}
		if (type[i] == 0) {
			std::cout << ">> Warning: Some data have undefined type!" << std::endl;
		}
	}

	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotRX[i] = tA;
 		t = tMin + i*tStep;

		hA = pPlot.AbsGEP(t);
		hB = pPlot.AbsGMP(t);
		plotRpY[i] = (hA/hB);

    }

	/// Graph

	TStyle *plain = new TStyle("Plain","Plain Style(no colors/fill areas)");
	plain->SetCanvasBorderMode(0);
	plain->SetPadBorderMode(0);
	plain->SetPadColor(0);
	plain->SetCanvasColor(0);
	plain->SetTitleColor(1);
	plain->SetStatColor(0);
	plain->SetTitleBorderSize(0);
	//gROOT->SetStyle(plain);
	plain->cd();
	Char_t const* title;
	Double_t x0, y0, s, w, h;
	x0 = 10;
	y0 = 10;
	s = 25;
	w = 900;
	h = 700;

	TCanvas* c = new TCanvas (uName("c",1), uName("Graph_",1), x0+s, y0+s, w, h);
	//c->SetLogy(); // logarithmic scale

	TMultiGraph *mgr1 = new TMultiGraph();

	TGraph *gr1 = new TGraph (nPoints, plotRX, plotRpY);
	gr1->SetLineColor(1);
	gr1->SetLineWidth(2);
	
	TGraphAsymmErrors *gr2 = new TGraphAsymmErrors (numE, dataX4, dataY4, 0, 0, dataU4, dataD4);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(.7);
	gr2->SetMarkerColor(2);

	TGraphAsymmErrors *gr3 = new TGraphAsymmErrors (numE, dataX4, dataY4, 0, 0, dataU4, dataD4);
	gr3->SetMarkerStyle(21);
	gr3->SetMarkerSize(.7);
	gr3->SetMarkerColor(4);

	mgr1->Add(gr1,"L");
	mgr1->Add(gr2,"P");
	mgr1->Add(gr3,"P");
	mgr1->Draw("A");

	mgr1->GetXaxis()->SetTitle("t [GeV^{2}]");
	mgr1->GetYaxis()->SetTitle("|G_E^p/G_M^p|");
	mgr1->GetXaxis()->CenterTitle();
	mgr1->GetYaxis()->CenterTitle();

	// Change the axis limits
	gPad->Modified();
	mgr1->GetXaxis()->SetLimits(3.,10.);
	mgr1->SetMinimum(0.01);
	mgr1->SetMaximum(3.);

	TLegend *lg1 = new TLegend (.5,.7,.87,.87);
	lg1->SetFillColor(kWhite);
 	lg1->AddEntry(gr1,title,"l");
	lg1->AddEntry(gr2,"World data with errors","ep");
	lg1->AddEntry(gr3,"World data with errors","ep");
	lg1->SetTextSize(0.03);
	lg1->Draw();
	
	/// End output

	theApp.Run(); //delete theApp;
	
	return EXIT_SUCCESS;
}

/**
 * @todo	split input data into arrays for plotting
 */