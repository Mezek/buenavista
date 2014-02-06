/**
 * $Date: 2013-06-18 10:03:42 +0200 (Tue, 18 Jun 2013) $
 * $Revision: 371 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/NucleonsH.cpp $
 * $Id: NucleonsH.cpp 371 2013-06-18 08:03:42Z bartos $
 *
 * @file
 * @brief	Main programe for nucleon form factors, v. hyperon.
 *  
 * <b>Compilation:</b>
 * @code
 * > g++ -o NucleonsH.exe NucleonsH.cpp `root-config --cflags` `root-config --libs` -lMinuit2
 * @endcode
 * 
 * <b>Usage</b> (to run as standalone application in console):
 * @code
 * > ./NucleonsH.exe --help
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
/** Global form factor data file.*/
char dataFile[] = "dataNucleonsWithMainzG.dat";     /**< Choose default data. */
char dataFile1[] = "dataNucleonsWithMainzGApprox.dat";
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

char parametersFile[] = "parNucleonsH.dat";                ///< Input parameters.
char outputFile[] = "outNucleonsH.dat";                    ///< Output parameters.
char xiFile[] = "outXiH.dat";                              ///< Chi2 output.
char matrixFile[] = "covarianceMatrixH.dat";               ///< Covariance matrix.
char debugFile[] = "outDebugH.dat";                        ///< Debug output.

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraphH.cpp"
#include "modules/HyperonsGetType.cpp"
#include "modules/Nucleon3GH.cpp"
#include "modules/NucleonFcnH.cpp"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserCovariance.h"
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

#include "NucleonsUsageH.cpp"
#include "NucleonsFitH.cpp"
#include "NucleonsChiH.cpp"
#include "NucleonsPlotH.cpp"
#include "NucleonsDebugH.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	time_t rawtime;
	time ( &rawtime );
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2013 Erik Bartoš" << std::endl;

	// Parse options

	performUsage(argc, argv);

	/// Start
	
	TApplication theApp("NucleonsH", &argc, argv);
	std::cout << "\n> Start  : " << ctime (&rawtime) << std::endl;
	
	/// Get data
	
	std::vector<double> param;
	std::vector<double> x, val, errUp, errDown;
	std::vector<int> type;

	/// Experimental points

	ExperimentalData A;
	//A.ReadData(globalArgs.data);
	A.ReadData(dataFile);
	std::cout << "\n>> World data:        `" << dataFile << "'" << std::endl;
	A.DataInfo();

	// Set precision
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

			// Fit: parameters, data

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
 * @todo	Check complex masses mwS2[i] for model.
 * @todo	Graphs.
 * @todo	Documentation.
 */
