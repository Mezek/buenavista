/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief       Main programe for pion form factors
 *  
 * <b>Compilation:</b>
 * @code
 * > g++ -o Pions.exe Pions.cpp `root-config --cflags` `root-config --libs` -lMinuit2
 * @endcode
 * 
 * <b>Usage</b> (to run as standalone application in console):
 * @code
 * > ./Pions.exe --help
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

#include "TROOT.h"
#include "TApplication.h"
#include "TComplex.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStopwatch.h"

char dataFile[] = "../Data/dataPi3Andrej.dat";         /**< Form factor data. */
//char dataFile0[] = "";
//char dataFile1[] = "";

//char parametersFile[] = "parNucleons-D.dat";
//char outputFile[] = "outNucleons-D.dat";
char parametersFile[] = "parPions.dat";                ///< Input parameters.
char outputFile[] = "outPions.dat";                    ///< Output parameters.
char xiFile[] = "outXi.dat";                           ///< Chi2 output.
char debugFile[] = "outDebug.dat";                     ///< Debug output.

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraph.cpp"
#include "modules/PionUam.cpp"
#include "modules/PionFcn.cpp"

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

#include "PionsUsage.cpp"
#include "PionsFit.cpp"
#include "PionsChi.cpp"
#include "PionsPlot.cpp"
#include "PionsDebug.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	time_t rawtime;
	time ( &rawtime );
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2014 Erik Bartoš" << std::endl;

	/// Parse options

	performUsage(argc, argv);

	/// Start
	
	TApplication theApp("Pions", &argc, argv);
	std::cout << "\n> Start  : " << ctime (&rawtime) << std::endl;
	
	/// Get data
	
	std::vector<double> param;
	std::vector<double> x, val, errUp, errDown;
	std::vector<int> type;

	/// Experimental points

	std::cout << "\n>> World data:" << std::endl;
	ExperimentalData A;
	A.ReadData(dataFile);	
	A.DataInfo();

	/// Set precision
	std::cout.precision(15); 
	//cout.setf(ios::scientific);

	switch (globalArgs.fit) {
		case 0: {

			/// Xi2: parameters, data

			performChi(globalArgs.parameters,globalArgs.data);

			/// Plot: parameters, data

			//performPlot(globalArgs.parameters,globalArgs.data);
			
			/// Debug: parameters, debugFile

			if ( globalArgs.verbose == 1) { performDebug(globalArgs.parameters,debugFile); }
			
			break;
		}
		case 1: {

			/// Fit: parameters, data

			performFit(globalArgs.parameters,globalArgs.data,globalArgs.output);

			performChi(globalArgs.output,globalArgs.data);

			//performPlot(globalArgs.output,globalArgs.data);

			if ( globalArgs.verbose == 1) { performDebug(globalArgs.output,debugFile); }

			break;
		}
		default: 
			break;
	}
	
	/// End output

	time ( &rawtime );
	std::cout << "\n> End    : " <<  ctime (&rawtime) << "=== " << std::endl;
	theApp.Run(); //delete theApp;
	
	return EXIT_SUCCESS;
}

/**
 * @todo	new fit, check chi squared
 * @todo	data, control at begining
 */
