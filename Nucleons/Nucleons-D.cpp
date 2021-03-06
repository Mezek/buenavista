/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief       Main programe for nucleon form factors, v. D.
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

#include "TROOT.h"
#include "TApplication.h"
#include "TComplex.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStopwatch.h"

///@{
char dataFile[] = "dataNucleonsDubnicka.dat";          /**< Form factor data. */
//char dataFile0[] = "dataNucleonsCleared.dat";
char dataFile1[] = "dataNucleonsApprox.dat";
char dataFile2[] = "dataNucleonsDubnicka.dat";
///@}

///@{
/** Single form factor data file.*/
char dataFileA[] = "../Data/dataProtonElectric.dat";
char dataFileB[] = "../Data/dataProtonMagnetic.dat";
char dataFileC[] = "../Data/dataNeutronElectric.dat";
char dataFileD[] = "../Data/dataNeutronMagnetic.dat";
char dataFileE[] = "../Data/dataProtonRatios.dat";
char dataFileF[] = "../Data/dataNeutronRatios.dat";
///@}

//char parametersFile[] = "parNucleons-D.dat";
//char outputFile[] = "outNucleons-D.dat";
char parametersFile[] = "parNucleons-Fit.dat";         ///< Input parameters.
char outputFile[] = "outNucleons-Fit.dat";             ///< Output parameters.
char xiFile[] = "outXi.dat";                           ///< Chi2 output.
char debugFile[] = "outDebug.dat";                     ///< Debug output.
char matrixFile[] = "covarianceMatrix.dat";            ///< Covariance matrix.
char radiusFile[] = "outRadius.dat";                   ///< Radius output.

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons-D.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraph.cpp"
#include "modules/Nucleon3G-D.cpp"
#include "modules/NucleonFcn.cpp"

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

#include "NucleonsUsage.cpp"
#include "NucleonsFit.cpp"
#include "NucleonsChi.cpp"
#include "NucleonsPlot.cpp"
#include "NucleonsDebug.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	time_t rawtime;
	time ( &rawtime );
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2013 Erik Bartoš" << std::endl;

	/// Parse options

	performUsage(argc, argv);

	/// Start
	
	TApplication theApp("Nucleons-D", &argc, argv);
	std::cout << "\n> Start  : " << ctime (&rawtime) << std::endl;
	
	/// Get data
	
	std::vector<double> param;
	std::vector<double> x, val, errUp, errDown;
	std::vector<int> type;

	/// Experimental points

	std::cout << "\n>> World data:" << std::endl;
	ExperimentalData A;
	A.ReadData(dataFile1);	
	A.ReadDataCrossSection(dataFile2);
	A.DataInfo();

	/// Set precision
	std::cout.precision(15); 
	//cout.setf(ios::scientific);

	switch (globalArgs.fit) {
		case 0: {

			/// Xi2: parameters, data

			performChi(globalArgs.parameters,globalArgs.data);

			/// Plot: parameters, data

			performPlot(globalArgs.parameters,globalArgs.data);
			
			/// Debug: parameters, debugFile

			if ( globalArgs.verbose == 1) { performDebug(globalArgs.parameters,debugFile); }
			
			break;
		}
		case 1: {

			/// Fit: parameters, data

			performFit(globalArgs.parameters,globalArgs.data,globalArgs.output);

			performChi(globalArgs.output,globalArgs.data);

			performPlot(globalArgs.output,globalArgs.data);

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
 * @todo	world data, control at begining
 */
