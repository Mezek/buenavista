/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief       Plot Eta.
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
#include "TStopwatch.h"

char dataFile1[] = "data/dataEtaPrime.dat";           ///< Form factor data.
char dataFile2[] = "data/dataEtaPrimeB.dat";
char parametersFile1[] = "parEtaPrimeFit.dat";        ///< Input parameters.
char parametersFile2[] = "parEtaPrimeB.dat";
char outputFile[] = "out-temp.dat";                   ///< Output parameters.

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraph.cpp"
#include "modules/MesonUam.cpp"
#include "modules/TransitionFcn.cpp"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserFcn.h"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2016 Erik Bartoš" << std::endl;

	/// Start
	
	TApplication theApp("EtaPrime", &argc, argv);

	/// Fit
	std::cout << "\n> EtaPrime fit:" << std::endl;
	std::cout << "> Input data:        `" << dataFile1 << "'" << std::endl;
	std::cout << "> Start parameters:  `" << parametersFile1 << "'" << std::endl;
	std::cout << "> Output parameters: `" << outputFile << "'" << std::endl;

	const int nPoints = 10000;
	double tMin;
	double tMax;
	double tStep;
	double tA;

	ExperimentalData Z;
	Z.ReadData(dataFile1);
	int nData = Z.size();
	std::cout << "\nNumber of fitted points:      " << nData << std::endl;
	std::vector<int> series = Z.Series();
	std::vector<std::string> names = Z.Name();
	std::vector<int> type = Z.Type();
	std::vector<double> x = Z.X();
	std::vector<double> val = Z.Val();
	std::vector<double> errUp = Z.ErrUp();
	std::vector<double> errDown = Z.ErrDown();
	std::vector<double> theta = Z.Theta();
	std::vector<double> energy = Z.Energy();

	int tParticle = 2;
	FFactorT trans(tParticle);
	int nPar = trans.numberOfParameters;
	trans.LoadParameters(parametersFile1);
	trans.PrintParameters();

	/// Create FCN function
	TransitionFcn fFCN(tParticle, Z.Type(), Z.X(), Z.Val(), Z.ErrUp(), Z.ErrDown(), Z.Theta(), Z.Energy());

	/// Minuit
	double step = 0.01;
	MnUserParameters upar;
	for (int i = 0; i < nPar; ++i) {
		upar.Add(trans.v[i].name, trans.v[i].val, step, trans.v[i].down, trans.v[i].up);
	}

	MnStrategy(2);
	/// Create minimizer
	MnMigrad migrad(fFCN, upar);
	MnMinimize minimize(fFCN, upar);
	MnSimplex simplex(fFCN, upar);
			
	//migrad.Fix(0.);
	//migrad.Fix(1.);
	//migrad.Fix(2.);
	//migrad.Fix(3.);
	//migrad.Fix(4.);
	//migrad.Fix(5.);

	/// Fitting
	TStopwatch timer;
	timer.Start();

	std::cout << "\n> Minimalization:" << std::endl;

	// operator(maxfcn,tolerance):
	// maxfcn: maximum function calls
	// migrad EDM: 0.001*tolerance*up
	// simplex:    tolerance*up
	
	std::cout << "> 1. minimalization..." << std::endl;	
	//FunctionMinimum min = migrad(100000,100.);

	//minimize.Fix(0.);
	//minimize.Fix(1.);
	//minimize.Fix(2.);
	//minimize.Fix(3.);
	//std::cout << "> ... 2. minimalization..." << std::endl;
	//FunctionMinimum min2 = migrad(10000,1.);

	std::cout << "> ... in progress..." << std::endl;
	//FunctionMinimum min9 = migrad(1000000,50.);
	FunctionMinimum min9 = minimize(1000000,50.);
	//FunctionMinimum min9 = simplex(1000000,1.);

	timer.Stop();

	/// Standard output
	std::cout << "> Minimum: " << min9 << std::endl;
	std::cout << "> FCN value:     " << min9.Fval() << std::endl;
	std::cout << "> FCN value/ndf: " << min9.Fval()/(nData-nPar) << std::endl;
	std::cout << "> Points       : " << nData << std::endl;

	std::cout << min9.UserState() << std::endl;
	std::cout << min9.UserCovariance() << std::endl;

	std::cout << "> Real time of minimalization: " << timer.RealTime() << " sec." << std::endl;

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
	std::ofstream os;
	os.open(outputFile);
	os.precision(8);
	for (int i = 0; i < nPar; ++i) {
		os << i+1 << ", " << min9.UserParameters().Name(i) << ", ";
		os.width(12);
		os << min9.UserParameters().Value(i) << ", ";
		os.width(12);
		os << min9.UserParameters().Error(i) << ", ";
		os << trans.v[i].down << ", " << trans.v[i].up << std::endl;
	}
	os.close();

	/// End output

	theApp.Run();
	//theApp.Terminate(0);
	
	return EXIT_SUCCESS;
}
