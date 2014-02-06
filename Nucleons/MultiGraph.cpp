/**
 * $Date: 2013-10-02 13:45:42 +0200 (Wed, 02 Oct 2013) $
 * $Revision: 383 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/MultiGraph.cpp $
 * $Id: MultiGraph.cpp 383 2013-10-02 11:45:42Z bartos $
 *
 * @file 
 * @brief	Plot single graph.
 *
 * <b>Compilation:</b>
 * @code
 * > g++ -o SingleGraph.exe SingleGraph.cpp -L/usr/lib64/root -I/usr/include/root `root-config --libs`.
 * @endcode
 * 
 * <b>Usage</b> (to run as standalone application in console):
 * @code
 * > ./SingleGraph.exe --help
 * @endcode
 */

#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <iomanip>
#include <time.h>

#include "TMath.h"
#include "TComplex.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TStopwatch.h"

char dataFile0[] = "../Data/dataNucleons.dat";
char dataFile1[] = "../Data/dataMainzG.dat";

#include "modules/ConstBasic.cpp"

//using namespace ROOT;

int main ( int argc, char **argv ) {

	TApplication theApp("Graph", &argc, argv);

	time_t rawtime;
	time ( &rawtime );
	TStopwatch timer;
	
	std::cout << "=== MultiGraph.cpp: version 1.0, (c) 2013 Erik Bartos ===" << std::endl;
	std::cout << "> Program: " << argv[0] << std::endl;
	std::cout << "> Start  : " << ctime (&rawtime) << std::endl;
	
	/// Get data
	std::vector<double> param;
	std::vector<double> x1, val1, errUp1, errDown1;
	std::vector<double> x2, val2, errUp2, errDown2;
	std::vector<int> type;
	char firstChar;
	std::string line;
	int num;

	/// Read source 1
	ifstream myDataFile (dataFile0);
	std::cout << ">> Data file: " << dataFile0 << std::endl;
	if (myDataFile.is_open()) {
		num = 0;
		double dA,dB,dC,dD;
		while (myDataFile.peek() != EOF) {
			firstChar = myDataFile.peek();
			if ( (firstChar == '%') || (firstChar == '#') || (firstChar == '*') || (firstChar == '/') ) {
				getline (myDataFile,line);
				//std::cout << line << std::endl;
				if ((firstChar == '#')) {
					std::cout << ">> " << line << std::endl;
				}
			}
			else {
				myDataFile >> dA >> dB >> dC >> dD;
				x1.push_back(dA);
				val1.push_back(dB);
				errUp1.push_back(dC);
				errDown1.push_back(dD);
				//std::cout << "> " << num << " " << x[num] << " " << val[num] << std::endl;
				getline (myDataFile,line);
				++num;
			}
		}
		myDataFile.close();
	}
	else std::cerr << ">> Error: Unable to open data file '" << dataFile0 << "'!" << std::endl;
	int nPoints1 = num;

	/// Read source 2
	ifstream myDataFile1 (dataFile1);
	std::cout << ">> Data file: " << dataFile1 << std::endl;
	if (myDataFile1.is_open()) {
		num = 0;
		double dA,dB,dC,dD;
		while (myDataFile1.peek() != EOF) {
			firstChar = myDataFile1.peek();
			if ( (firstChar == '%') || (firstChar == '#') || (firstChar == '*') || (firstChar == '/') ) {
				getline (myDataFile1,line);
				//std::cout << line << std::endl;
				if ((firstChar == '#')) {
					std::cout << ">> " << line << std::endl;
				}
			}
			else {
				myDataFile1 >> dA >> dB >> dC >> dD;
				x2.push_back(dA);
				val2.push_back(dB);
				errUp2.push_back(dC);
				errDown2.push_back(dD);
				//std::cout << "> " << num << " " << x[num] << " " << val[num] << std::endl;
				getline (myDataFile1,line);
				++num;
			}
		}
		myDataFile1.close();
	}
	else std::cerr << ">> Error: Unable to open data file '" << dataFile0 << "'!" << std::endl;
	int nPoints2 = num;
	
	/// Set precision
	std::cout.precision(20); 
	//cout.setf(ios::scientific); 

	/// Plot graph
	double plotX1[nPoints1], plotY1[nPoints1];
	double plotEX1[nPoints1], plotEY1[nPoints1];
	for (int i=0; i < nPoints1; ++i) {
		plotX1[i] = x1[i];
		plotEX1[i] = x1[i];
		plotY1[i] = val1[i];
		plotEY1[i] = errUp1[i];
	}

	double plotX2[nPoints2], plotY2[nPoints2];
	double plotEX2[nPoints2], plotEY2[nPoints2];
	for (int i=0; i < nPoints2; ++i) {
		plotX2[i] = x2[i];
		plotEX2[i] = x2[i];
		plotY2[i] = val2[i];
		plotEY2[i] = errUp2[i];
	}

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


	TCanvas *m1 = new TCanvas ("m1", "Graph", 10, 10, 700, 500);
    m1->SetLogy(); // logarithmic scale

	TMultiGraph *mgr1 = new TMultiGraph();

    TGraph *gr1 = new TGraph (nPoints1, plotX1, plotY1);
	//gr1->Draw("AP");
	gr1->SetMarkerStyle(21);
	gr1->SetMarkerSize(.3);
	gr1->SetTitle("#sigma_{exp}");

	TGraph *gr2 = new TGraph (nPoints2, plotX2, plotY2);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(.7);
	gr2->SetMarkerColor(2);

	mgr1->Add(gr1,"P");
	mgr1->Add(gr2,"P");
	mgr1->Draw("A");

	mgr1->GetXaxis()->SetTitle("t [GeV^{2}]");
	mgr1->GetYaxis()->SetTitle("|G|");
    mgr1->GetXaxis()->CenterTitle();
    mgr1->GetYaxis()->CenterTitle();

	// Change the axis limits
	gPad->Modified();
	mgr1->GetXaxis()->SetLimits(-7.5,15.);
	mgr1->SetMinimum(0.001);
	mgr1->SetMaximum(10.);

	TLegend *lg1 = new TLegend (.5,.7,.87,.87);
	lg1->SetFillColor(kWhite);
 	lg1->AddEntry(gr1,"|G|","p");
	lg1->AddEntry(gr2,"World data with errors","ep");
	lg1->SetTextSize(0.03);
	lg1->Draw();
	
	time ( &rawtime );
	std::cout << "\n> End of program: " <<  ctime (&rawtime) << std::endl;
	theApp.Run();
	
	return 0;
}

/**
 * @todo	...
 */
