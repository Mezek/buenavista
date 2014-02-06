/**
 * $Date: 2013-10-02 13:45:42 +0200 (Wed, 02 Oct 2013) $
 * $Revision: 383 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/Nucleons.cpp $
 * $Id: Nucleons.cpp 383 2013-10-02 11:45:42Z bartos $
 *
 * @file
 * @brief       Main programe for nucleon form factors.
 *  
 * <b>Compilation:</b>
 * @code
 * > g++ -o Nucleons.exe Nucleons.cpp `root-config --cflags` `root-config --libs` -lMinuit2
 * @endcode
 * 
 * <b>Usage</b> (to run as standalone application in console):
 * @code
 * > ./Nucleons.exe --help
 * @endcode
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

char dataFile0[] = "dataNucleons.dat";
//char dataFile1[] = "dataNucleonsCleared.dat";
char dataFile1[] = "dataNucleonsApprox.dat";
char dataFileA[] = "../Data/dataProtonElectric.dat";
char dataFileB[] = "../Data/dataProtonMagnetic.dat";
char dataFileC[] = "../Data/dataNeutronElectric.dat";
char dataFileD[] = "../Data/dataNeutronMagnetic.dat";
char dataFileE[] = "../Data/dataProtonRatios.dat";
char dataFileF[] = "../Data/dataNeutronRatios.dat";

char parametersFile[] = "parNucleons-Fit.dat";
char outputFile[] = "outNucleons-Fit.dat";
char xiFile[] = "outXi.dat";

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons-D.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraph.cpp"
#include "modules/Nucleon3G-D.cpp"
#include "modules/NucleonFcn.cpp"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnPrint.h"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	TApplication theApp("Nucleons", &argc, argv);

	time_t rawtime;
	time ( &rawtime );
	TStopwatch timer;
	
	std::cout << "=== Nucleons.cpp: version 1.0, (c) 2013 Erik Bartos ===" << std::endl;
	std::cout << "> Program: " << argv[0] << std::endl;
	std::cout << "> Start  : " << ctime (&rawtime) << std::endl;
	
	/// Get data
	int nData = 0;

	/*std::vector<double> param;
	std::vector<double> x, val, errUp, errDown;
	std::vector<int> type;*/

	/// Experimental points
	ExperimentalData A,B,C,D,E,F,W,Z;

	A.ReadData(dataFileA);
	B.ReadData(dataFileB);
	C.ReadData(dataFileC);
	D.ReadData(dataFileD);
	E.ReadData(dataFileE);
	F.ReadData(dataFileF);
	
	int numA = A.size();
	int numB = B.size();
	int numC = C.size();
	int numD = D.size();
	int numE = E.size();
	int numF = F.size();
	int numAll = numA+numB+numC+numD+numE+numF;
	
	std::cout << "Data protonElectric:           " << numA << std::endl;
	std::cout << "Data protonMagnetic:           " << numB << std::endl;
	std::cout << "Data neutronElectric:          " << numC << std::endl;
	std::cout << "Data neutronMagnetic:          " << numD << std::endl;
	std::cout << "Data protonRatios:             " << numE << std::endl;
	std::cout << "Data neutronRatios:            " << numF << std::endl;
	std::cout << "Number of experimental points: " << numAll << std::endl;

	Z.ReadData(dataFile1);
	nData = Z.size();
	std::cout << "\nNumber of used points:         " << nData << "\t(" << 100.*nData/numAll << "%)" << std::endl;
	Z.CheckData();

	W.ReadData(dataFile0);
	std::vector<double> x = W.X();
	std::vector<double> val = W.Val();
	std::vector<double> errUp = W.ErrUp();
	std::vector<double> errDown = W.ErrDown();
	std::vector<int> type = W.Type();
		
	/// Set precision
	std::cout.precision(20); 
	//cout.setf(ios::scientific); 

	const int fit = 0;
	
	switch (fit) {
		case 0:
			break;
		case 1: {
	
	/// Fit of form factor
	FFactor pFF(12);
	//~pFF.LoadParameters(pFile);
	//~std::cout << "> Number of parameters: " << pFF.numberOfParameters << std::endl;	
	//~pFF.PrintParameters();

	pFF.LoadParameters(parametersFile);
	int nPar = pFF.numberOfParameters;
	std::cout << "> Number of FF parameters: " << nPar << std::endl;
	//pFF.PrintParameters();

	/// Create FCN function
	NucleonFcn fFCN(nPar, Z.Type(), Z.X(), Z.Val(), Z.ErrUp(), Z.ErrDown(), Z.Theta(), Z.Energy());

	// Minuit
	double step = 0.01;
	MnUserParameters upar;
	upar.Add("t_1s  ", pFF.v[0].val, step, 0.9, 5.0);
	upar.Add("t_1v  ", pFF.v[1].val, step, 2.0, 2.5);
	upar.Add("t_2s  ", pFF.v[2].val, step, 1.0, 2.0);
	upar.Add("t_2v  ", pFF.v[3].val, step, 1.8, 2.5);
	upar.Add("f_Om  ", pFF.v[4].val, step, -1.7, .1);
	upar.Add("f_Ph  ", pFF.v[5].val, step, -0.7, 2.8);
	upar.Add("f_Om1 ", pFF.v[6].val, step, -1, 1.);
	upar.Add("f_Ph1 ", pFF.v[7].val, step, -0.9, 0.);
	upar.Add("f_Rh  ", pFF.v[8].val, step, -2.7, 0.);
	upar.Add("f_Om_T", pFF.v[9].val, step, -2.7, 1.7);
	upar.Add("f_Ph_T", pFF.v[10].val, step, -3.2, -1.0);
	upar.Add("f_Om1T", pFF.v[11].val, step, 0., 1.3);

	/// Create Migrad minimizer
	MnMigrad migrad(fFCN, upar, 1);

	/*migrad.Fix(0.);
	migrad.Fix(1.);
	migrad.Fix(2.);
	migrad.Fix(3.);
	/*migrad.Fix(4.);
	migrad.Fix(5.);
	migrad.Fix(6.);
	migrad.Fix(7.);
	migrad.Fix(8.);
	migrad.Fix(9.);
	migrad.Fix(10.);
	migrad.Fix(11.);
	migrad.Fix(12.);
	migrad.Fix(13.);
	migrad.Fix(14.);
	migrad.Fix(15.);*/

	/// Fitting
	timer.Start();

	//std::cout << "> 1. minimalization..." << std::endl;	
	//FunctionMinimum min = migrad(100000,1.);

	//.migrad.Fix(15.);
	//std::cout << "> ... 2. minimalization..." << std::endl;
	//FunctionMinimum min2 = migrad(10000,1.);
	
	std::cout << "> ... last minimalization..." << std::endl;
	FunctionMinimum min9 = migrad(100000,1.);

	timer.Stop();

	/// Standard output
	std::cout << "> Minimum: " << min9 << std::endl;
	std::cout << "> FCN value: " << min9.Fval() << std::endl;
	std::cout << "> Real time of minimalization: " << timer.RealTime() << " sec." << std::endl;

	//pFF.PrintParameters();
	
	//~ // Scan and Plot parameters
 	//~ {
    //~ MnScan scan(fFCN, upar, 1);
    //~ std::cout << "Scan parameters: " << scan.Parameters() << std::endl;
    //~ MnPlot plot;
    //~ for (unsigned int i = 0; i < upar.VariableParameters(); ++i) {
		//~ std::vector<std::pair<double, double> > xy = scan.Scan(i);
		//~ //std::vector<std::pair<double, double> > xy = scan.scan(0);
		//~ std::cout << i+1 << ". parameter - `" << upar.Name(i) << "'"<< std::endl;
		//~ plot(xy);
	//~ }
    //~ std::cout << scan.Parameters() << std::endl;
	//~ }

	/// Output parameters
	ofstream myOutputParam (outputFile);
	for (int i = 0; i < nPar; ++i) {
		myOutputParam << min9.UserParameters().Value(i) << std::endl;
		//std::cout << min9.UserParameters().Value(i) << std::endl;
	}
	myOutputParam.close();

		break;
	}
		default: {}
	}

	/*
	// Output parcial Xi2
	ofstream myOutputXi (xiFile);
	FFactor bornTBW(12);
	bornTBW.LoadParameters(outputFile);
	double chi = 0.;
	double chi2 = 0.;
	double Chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;
	TComplex nFv,nE,nM;
	for(int n = 0; n < Z.size(); ++n) {
		// "protonElectric"
		if ((Z.type[n] == 1) || (Z.type[n] == 2)) {
			nFv = nFv.Abs(bornTBW.GEP(Z.x[n]));
		}
		// "protonMagnetic"
		if ((Z.type[n] == 3) || (Z.type[n] == 4)) {
			nFv = nFv.Abs(bornTBW.GMP(Z.x[n]));
		}
		// "neutronElectric"
		if ((Z.type[n] == 5) || (Z.type[n] == 6)) {
			nFv = nFv.Abs(bornTBW.GEN(Z.x[n]));
		}
		// "neutronMagnetic"
		if ((Z.type[n] == 7) || (Z.type[n] == 8)) {
			nFv = nFv.Abs(bornTBW.GMN(Z.x[n]));
		}
		// "protonRatios"
		if (Z.type[n] == 9) {
			nE = bornTBW.GEP(Z.x[n]);
			nM = bornTBW.GMP(Z.x[n]);
			nFv = (1.+ammP)*nFv.Abs(nE/nM);
		}
		// "neutronRatios"
		if (Z.type[n] == 10) {
			nE = bornTBW.GEN(Z.x[n]);
			nM = bornTBW.GMN(Z.x[n]);
			nFv = (ammN)*nFv.Abs(nE/nM);
		}

		delta = nFv.Re()-Z.val[n];
		sigmaAverage = (Z.errUp[n] + Z.errDown[n])/2.;
		chi = delta/sigmaAverage;
		chi2 = chi*chi;
		Chi2 += chi2;
		if ( chi2 > 7.) {
			myOutputXi << Z.type[n] << "\t" << chi2 << "\t" << n << ". " << Z.x[n] << " " << Z.val[n] << std::endl;
		}
	}
	myOutputXi << "Chi2: " << Chi2 << std::endl;
	myOutputXi.close();
	std::cout << "\n> Chi2: " << Chi2 << std::endl;
	*/
	
	/// Plot FFs with data points
	FFactor nPlot(12);
	nPlot.LoadParameters(outputFile);
	//nPlot.PrintParameters();
	nPlot.CheckFormFactor("all", 0.01);

	const int nPoints = 2500;
	const int nC = 500;
	double plotX[nPoints], plotY0[nPoints], plotY1[nPoints], plotY2[nPoints], plotY3[nPoints];
	double dataX0[numA], dataX1[numB], dataX2[numC], dataX3[numD];
	double dataY0[numA], dataY1[numB], dataY2[numC], dataY3[numD];
	//std::fill(dataX3, dataX3+num3, 0);
	double tMin = -4.5;
	double tMax = 15.0;
	double tStep = (tMax-tMin)/nPoints;
	double tA = tMin;
	double res;
	TComplex z0(tMin,0.000001);
	TComplex z,zY,zW,gA,gB,gC,gD;

	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotX[i] = tA;
		//~???
 		//~if (tA < t0v) { 
 		z = tMin + i*tStep;
 		//~}
 		//~else {
 			//~z = z0 + i*shag;
 		//~}
		gA = nPlot.ScalarOne(z);
		gB = nPlot.VectorOne(z);
		gC = nPlot.ScalarTwo(z);
		gD = nPlot.VectorTwo(z);

		plotY0[i] = TComplex::Abs(gA+gB + z/(4*massP*massP)*(gC+gD));
		plotY1[i] = TComplex::Abs(gA+gB + gC+gD);
		plotY2[i] = TComplex::Abs(gA-gB + z/(4*massN*massN)*(gC-gD));
		plotY3[i] = TComplex::Abs(gA-gB + gC-gD);
    }
	
	int k0 = 0 , k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	for (int i = 0; i < W.size(); ++i) {
		// "protonElectric"
		if ((type[i] == 1) || (type[i] == 2)) {
			dataX0[k0] = x[i];
			dataY0[k0] = val[i];
			//std::cout << i << ": " << x[i] << " " << val[i] << std::endl;
			//std::cout << i << ": " << dataX0[k0] << " " << dataY0[k0] << std::endl;
			++k0;
		}
		// "protonMagnetic"
		if ((type[i] == 3) || (type[i] == 4)) {
			dataX1[k1] = x[i];
			dataY1[k1] = val[i];
			++k1;
		}
		// "neutronElectric"
		if ((type[i] == 5) || (type[i] == 6)) {
			dataX2[k2] = x[i];
			dataY2[k2] = val[i];
			++k2;
		}
		// "neutronMagnetic"
		if ((type[i] == 7) || (type[i] == 8)) {
			dataX3[k3] = x[i];
			dataY3[k3] = val[i];
			//std::cout << " " << k3 << ". " << dataX3[k3] << std::endl;
			++k3;
		}
		if ((type[i] == 9) || (type[i] == 10)) {
			//std::cout << " " << k4 << ". ";
			++k4;
		}
		if (type[i] == 11) {
			std::cout << ">> Warning: No type of form factor defined!" << std::endl;
		}
	}
	
	
	PlotGraph graf1(10);
	const Char_t* title = "|G_{M}^{n}| with data";
	//graf1.viewPlusData(nPoints,plotX,plotY0,k0,dataX0,dataY0,title);
	//graf1.viewPlusData(nPoints,plotX,plotY1,k1,dataX1,dataY1,title);
	//graf1.viewPlusData(nPoints,plotX,plotY2,k2,dataX2,dataY2,title);
	graf1.viewPlusData(nPoints,plotX,plotY3,k3,dataX3,dataY3,title);
	
	//PlotGraph graf2;
	//graf2.view4(nPoints,plotX,plotY0,plotY1,plotY2,plotY3);
	//graf2.view4Exp(nPoints,plotX,plotY0,plotY1,plotY2,plotY3,k0,dataX0,dataY0,k1,dataX1,dataY1,k2,dataX2,dataY2,k3,dataX3,dataY3);

	/// Plot ratios
	FFactor pPlot(12);
	pPlot.LoadParameters(outputFile);

	double plotRX[nPoints], plotRY[nPoints];
	double dataX4[numE], dataX5[numF];
	double dataY4[numE], dataY5[numF];
	double dataEU4[numE], dataED4[numE];
	tMin = 0.0;
	tMax = 19.0;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	TComplex w,hA,hB;

	for (int i = 0; i < numE; ++i) {
		dataX4[i] = -E.X()[i];
		dataY4[i] = E.Val()[i];
		dataEU4[i] = 0.;
		dataED4[i] = E.ErrUp()[i];
		//std::cout << i << ". " << dataX4[i] << " " << dataY4[i] <<  std::endl;
	}
	
	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotRX[i] = tA;
 		w = tMin + i*tStep;

		hA = pPlot.GEP(-w);
		hB = pPlot.GMP(-w);

		plotRY[i] = (1+ammP)*hA/hB;
		//std::cout << i << ". " << w << " " << hA <<  std::endl;
    }

	PlotGraph grafR(13);
	const Char_t* titleR = "mu_p*G_E^p/G_M^p";
	grafR.viewPlusDataE(nPoints,plotRX,plotRY,numE,dataX4,dataY4,dataEU4,dataED4,titleR);


	/*/// Plot isoFF
	FFactor oPlot(12);
	oPlot.LoadParameters(outputFile);
	oPlot.PerformCheck("all", 0.001);
	
	double plotIX[nPoints], plotIY[nPoints];
	tMin = -0.5;
	tMax = 0.0;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	TComplex hA;
	
	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotIX[i] = tA;
 		z = tMin + i*tStep;

		//hA = oPlot.ScalarOne(z);
		//hA = nPlot.VectorOne(z);
		hA = nPlot.ScalarTwo(z);
		//hA = nPlot.VectorTwo(z);

		plotIY[i] = TComplex::Abs(hA);
		//std::cout << i << ". " << z << " " << hA <<  std::endl;
    }

	PlotGraph grafI;
	const Char_t* tit = "IsoFF";
	grafI.view(nPoints,plotIX,plotIY,tit);
	*/

	/// End output

	time ( &rawtime );
	std::cout << "\n> End of program: " <<  ctime (&rawtime) << std::endl;
	theApp.Run(); //delete theApp;
	
	return EXIT_SUCCESS;
}

/**
 * @todo	print GM_p and compare with D. graph
 * @todo	new fit, check chi squared
 * @todo	perform D. fit with fortran minuit
 */
