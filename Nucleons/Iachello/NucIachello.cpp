/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Main programe for nucleon form factors, Iachello model, v. with cross section data.
 *  
 * <b>Compilation:</b>
 * @code
 * > g++ -o NucIachello.exe NucIachello.cpp `root-config --cflags` `root-config --libs` -lMinuit2
 * @endcode
 * 
 * <b>Usage</b> (to run as standalone application in console):
 * @code
 * > ./NucIachello.exe --help
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
//char dataFile[] = "../../Data/dataProtonRatios.dat";
char dataFile1[] = "../../Data/dataProtonRatios.dat";  /**< Choose default data. */
char dataFile2[] = "../../Data/dataEmpty.dat";
///@}

///@{
/** Single form factor data file.*/
char dataFileA[] = "../../Data/dataProtonElectric.dat";
char dataFileB[] = "../../Data/dataProtonMagnetic.dat";
char dataFileC[] = "../../Data/dataNeutronElectric.dat";
char dataFileD[] = "../../Data/dataNeutronMagnetic.dat";
char dataFileE[] = "../../Data/dataProtonRatios.dat";
char dataFileF[] = "../../Data/dataNeutronRatios.dat";
///@}

char parametersFile[] = "parNucIachello.dat";          ///< Input parameters.
char outputFile[] = "outNucIachello.dat";              ///< Output parameters.
char xiFile[] = "outXi.dat";                           ///< Chi2 output.
char matrixFile[] = "covarianceMatrix.dat";            ///< Covariance matrix.
char radiusFile[] = "outRadius.dat";                   ///< Radius output.
char debugFile[] = "outDebug.dat";                     ///< Debug output.
char tableFile[] = "outTable.dat";                     ///< Tabulation. 

#include "../modules/ConstBasic.cpp"
#include "../modules/ConstMesons.cpp"
#include "../modules/ExperimentalData.cpp"
#include "../modules/NucleonIachello.cpp"
#include "../modules/NucleonFcnCS.cpp"
#include "PlotGraphIachello.cpp"

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

#include "NucIachelloUsage.cpp"
#include "NucIachelloChi.cpp"
#include "NucIachelloPlot.cpp"
#include "NucIachelloFit.cpp"
#include "NucIachelloDebug.cpp"
#include "NucIachelloTable.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	time_t rawtime;
	time ( &rawtime );
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2013 Erik Bartoš" << std::endl;

	/// Parse options

	performUsage(argc, argv);

	/// Start
	
	TApplication theApp("NucleonsIachello", &argc, argv);
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

			performChi(globalArgs.parameters,globalArgs.ffdata,globalArgs.csdata);

			/// Debug: parameters, debugFile

			if ( globalArgs.verbose == 1) {
				//performRadius(globalArgs.parameters,matrixFile);
				performDebug(globalArgs.parameters,debugFile);
			}

			/// Plot: parameters, data

			//performPlot(globalArgs.parameters,globalArgs.ffdata,globalArgs.csdata);

			if ( globalArgs.tabulate == 1) { performTable(globalArgs.parameters,tableFile); }
				
			break;
		}
		case 1: {

			/// Fit: parameters, data

			performFit(globalArgs.parameters,globalArgs.ffdata,globalArgs.csdata,globalArgs.output);

			performChi(globalArgs.output,globalArgs.ffdata,globalArgs.csdata);

			performPlot(globalArgs.output,globalArgs.ffdata,globalArgs.csdata);

			if ( globalArgs.verbose == 1) { performDebug(globalArgs.output,debugFile); }
			if ( globalArgs.tabulate == 1) { performTable(globalArgs.output,tableFile); }

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
 * @todo       Graphs.
 * @todo       Documentation.
 */
