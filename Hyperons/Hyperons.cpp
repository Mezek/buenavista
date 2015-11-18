/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief       Main programe for hyperon form factors.
 *  
 * <b>Compilation:</b>
 * @code
 * > g++ -o Hyperons.exe Hyperons.cpp `root-config --cflags` `root-config --libs` -lMinuit2
 * @endcode
 * 
 * <b>Usage</b> (to run as standalone application in console):
 * @code
 * > ./Hyperons.exe --help
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
#include "TStopwatch.h"

///@{
/** Global form factor data file.*/
char dataFile[] = "dataNucleonsCleared.dat";
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

char parametersFile[] = "parHyperons.dat";             ///< Input parameters.
char outputFile[] = "outHyperons.dat";                 ///< Output parameters.
char xiFile[] = "outXiHyperons.dat";                   ///< Chi2 output.
char debugFile[] = "outDebugHyperons.dat";             ///< Debug output.

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraphH.cpp"
#include "modules/HyperonsGetType.cpp"
#include "modules/Hyperon3G.cpp"
#include "modules/NucleonFcnH.cpp"

#include "HyperonsUsage.cpp"
#include "HyperonsPlot.cpp"
#include "HyperonsRadii.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	TApplication theApp("Hyperons", &argc, argv);

	time_t rawtime;
	time ( &rawtime );
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2013 Erik Bartoš" << std::endl;

	/// Parse options

	performUsage(argc, argv);
	
	/// Start
	
	std::cout << "\n> Start  : " << ctime (&rawtime) << std::endl;

	std::vector<double> paramHyp;

	// Set precision
	std::cout.precision(20); 
	//cout.setf(ios::scientific);

	performPlot(globalArgs.parameters,globalArgs.data);

	//performRadii(globalArgs.parameters,globalArgs.output);

	/// Debug: parameters, debugFile

	//if ( globalArgs.verbose == 1) { performDebug(globalArgs.parameters,debugFile); }

	/// End output

	time ( &rawtime );
	std::cout << "\n> End    : " <<  ctime (&rawtime) << "=== " << std::endl;
	theApp.Run();
	
	return EXIT_SUCCESS;
}
