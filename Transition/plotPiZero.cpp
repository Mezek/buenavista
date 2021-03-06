/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief       Plot PiZero.
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

char dataFile1[] = "data/dataPiZero.dat";            ///< Form factor data.
char parametersFile1[] = "parPiZeroFit.dat";         ///< Input parameters.
char parametersFile2[] = "parPiZeroFit.dat";
char outputFile[] = "outEta-temp.dat";               ///< Output parameters.

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraph.cpp"
#include "modules/MesonUam.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2016 Erik Bartoš" << std::endl;

	/// Start
	
	TApplication theApp("FF", &argc, argv);

	/// Plot
	std::cout << "\n> Plotting:" << std::endl;
	std::cout << "> Plotted parameters:          `" << parametersFile1 << "'" << std::endl;
	std::cout << "> Plotted form factor data:    `" << dataFile1 << "'" << std::endl;

	const int nPoints = 10000;
	double tMin;
	double tMax;
	double tStep;
	double tA;

	ExperimentalData W;
	W.ReadData(dataFile1);
	std::cout << "\nNumber of plotted points:      " << W.size() << std::endl;
	std::vector<int> series = W.Series();
	std::vector<std::string> names = W.Name();
	std::vector<int> type = W.Type();
	std::vector<double> x = W.X();
	std::vector<double> val = W.Val();
	std::vector<double> errUp = W.ErrUp();
	std::vector<double> errDown = W.ErrDown();
	std::vector<double> theta = W.Theta();
	std::vector<double> energy = W.Energy();

	// Plot function
	
	FFactorT trans(3);
	trans.LoadParameters(parametersFile1);
	trans.PrintParameters();

	// Debug
	std::cout << "\n>> Debug: " << std::endl;
	double t1 = trans.A(0);
	double t2 = trans.A(1);
	double et1 = trans.E(0);
	double et2 = trans.E(1);

	double h1 = sqrt(t1/t0t - 1.);
	double eh1 = et1/(2.*t1*h1);
	double h2 = sqrt(t2/t0t - 1.);
	double eh2 = et2/(2.*t2*h2);
	std::cout << ">> q_in^S : "<< h1 << " +/- " << eh1 << std::endl;
	std::cout << ">> q_in^V : "<< h2 << " +/- " << eh2 << std::endl;

	double plotRX[nPoints], plotRpY[nPoints], plotRSY[nPoints], plotRVY[nPoints];
	
	tMin = -13.0;
	tMax = 5.0;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	TComplex t,hA,hB,hC,hD;

	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotRX[i] = tA;
 		t = tMin + i*tStep;

		hA = trans.FFAbsVal(t);
		plotRpY[i] = hA;
		plotRSY[i] = trans.ScalarP(t).Rho();
		plotRVY[i] = trans.VectorP(t).Rho();
    }

	// Graph

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
	Double_t x0, y0, s, w, h;
	x0 = 10;
	y0 = 10;
	s = 25;
	w = 900;
	h = 700;

	TCanvas* c = new TCanvas (uName("c",1), uName("Graph_",1), x0+s, y0+s, w, h);
	c->SetLogy(); // logarithmic scale

	TMultiGraph *mgr1 = new TMultiGraph();

	TGraph *gr1 = new TGraph (nPoints, plotRX, plotRpY);
	gr1->SetLineColor(1);
	gr1->SetLineWidth(2);
	mgr1->Add(gr1,"L");

	TGraph *gr2 = new TGraph (nPoints, plotRX, plotRSY);
	gr2->SetLineColor(2);
	gr2->SetLineWidth(2);
	//mgr1->Add(gr2,"L");

	TGraph *gr3 = new TGraph (nPoints, plotRX, plotRVY);
	gr3->SetLineColor(4);
	gr3->SetLineWidth(2);
	//mgr1->Add(gr3,"L");
	
	int numS = series.size();
	TGraphAsymmErrors *g[numS];
	
	for (int i = 0; i < numS; i++) {
		std::vector<double> dataX;
		std::vector<double> dataY;
		std::vector<double> dataU;
		std::vector<double> dataD;

		for (int j = 0; j < W.size(); j++) {
			if (series[i] == type[j]) {
				dataX.push_back(x[j]);
				dataY.push_back(val[j]);
				dataU.push_back(errUp[j]);
				dataD.push_back(errDown[j]);
			}
			if (type[j] == 0) {
				std::cout << ">> Warning: Some data have undefined type!" << std::endl;
			}
		}
		g[i] = new TGraphAsymmErrors (dataX.size(), &(dataX[0]), &(dataY[0]), 0, 0, &(dataU[0]), &(dataD[0]));
		mgr1->Add(g[i],"P");

		g[i]->SetMarkerStyle(20+i);
		g[i]->SetMarkerSize(1);
		g[i]->SetMarkerColor(2+i);
	}
	mgr1->Draw("A");

	mgr1->GetXaxis()->SetTitle("t [GeV^{2}]");
	mgr1->GetYaxis()->SetTitle("|F_{#pi^{0}}|");
	mgr1->GetXaxis()->CenterTitle();
	mgr1->GetYaxis()->CenterTitle();

	// Change the axis limits
	gPad->Modified();
	mgr1->GetXaxis()->SetLimits(-14.,5.);
	//mgr1->SetMinimum(0.01);
	//mgr1->SetMaximum(3.);

	Char_t const* title = "U&A model";

	TLegend *lg1 = new TLegend (.17,.67,.47,.87);
	lg1->SetFillColor(kWhite);
 	lg1->AddEntry(gr1,title,"l");
	//lg1->AddEntry(gr2,"Kelly","l");
	for (int i = 0; i < numS; i++) {
		lg1->AddEntry(g[i],names[i].c_str(),"ep");
	}
	lg1->SetTextSize(0.03);
	lg1->Draw();

	c->Update();

	c->SaveAs("gPiZero.pdf");

	/// End output

	theApp.Run();
	//theApp.Terminate(0);
	
	return EXIT_SUCCESS;
}
