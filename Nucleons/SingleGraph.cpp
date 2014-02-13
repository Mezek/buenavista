/**
 * $Date$
 * $Revision$
 * $Author$
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

//char dataFile0[] = "../Data/dataMainz.dat";
char dataFile0[] = "../Data/dataMainzG.dat";

#include "modules/ConstBasic.cpp"

//using namespace ROOT;

int main ( int argc, char **argv ) {

	TApplication theApp("Graph", &argc, argv);

	time_t rawtime;
	time ( &rawtime );
	TStopwatch timer;
	
	std::cout << "=== SingleGraph.cpp: version 1.0, (c) 2013 Erik Bartos ===" << std::endl;
	std::cout << "> Program: " << argv[0] << std::endl;
	std::cout << "> Start  : " << ctime (&rawtime) << std::endl;
	
	/// Get data
	std::vector<double> param;
	std::vector<double> x, val, errUp, errDown;
	std::vector<int> type;
	char firstChar;
	std::string line;
	int num;

	/// Read source
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
				x.push_back(dA);
				val.push_back(dB);
				errUp.push_back(dC);
				errDown.push_back(dD);
				//std::cout << "> " << num << " " << x[num] << " " << val[num] << std::endl;
				getline (myDataFile,line);
				++num;
			}
		}
		myDataFile.close();
	}
	else std::cerr << ">> Error: Unable to open data file '" << dataFile0 << "'!" << std::endl;
	int nPoints = num;
		
	/// Set precision
	std::cout.precision(20); 
	//cout.setf(ios::scientific); 

	/// Plot graph

	double plotX[nPoints], plotY[nPoints];
	for (int i=0; i < nPoints; ++i) {
		plotX[i] = x[i];
		plotY[i] = val[i];
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
		
	TCanvas *s1 = new TCanvas ("s1", "Graph", 10, 10, 700, 500);
    //s1->SetLogy(); // logarithmic scale
	
	TGraph *gr1 = new TGraph (nPoints, plotX, plotY);
	gr1->Draw("AP");
	gr1->SetMarkerStyle(21);
	gr1->SetMarkerSize(.3);
	gr1->SetTitle("#sigma_{exp}");
	
	time ( &rawtime );
	std::cout << "\n> End of program: " <<  ctime (&rawtime) << std::endl;
	theApp.Run();
	
	return 0;
}

/**
 * @todo	...
 */
