/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief       Main programe for nucleon form factors, v. MC.
 *  
 * <b>Compilation:</b>
 * @code
 * > g++ -o NucleonsMC.exe NucleonsMC.cpp `root-config --cflags` `root-config --libs` -lMinuit2.
 * @endcode
 * 
 * <b>Usage</b> (to run as standalone application in console):
 * @code
 * > ./NucleonsMC.exe --help
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

// 'A' = protonElectric
// 'B' = protonMagnetic
// 'C' = neutronElectric
// 'D' = neutronMagnetic

char dataFile1[] = "dataNucleonsApprox.dat";
char dataFile2[] = "dataMainzRatioFit.dat";

char parametersFile[] = "parNucleons-Mainz.dat";
char outputFile[] = "outNucleons-MC.dat";
char radiusFile[] = "outRadius-MC.dat";

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons-D.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraph.cpp"
#include "modules/Nucleon3G-D.cpp"
#include "modules/NucleonFcnCS.cpp"

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

#include "NucleonsUsageMC.cpp"
#include "NucleonsFitMC.cpp"
#include "NucleonsRadiusMC.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	time_t rawtime;
	time ( &rawtime );
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2013 Erik Bartoš" << std::endl;

	/// Parse options

	performUsage(argc, argv);	

	/// Start
	
	std::cout << "\n> Start  : " << ctime (&rawtime) << std::endl;
	
	/// Set precision
	std::cout.precision(15); 
	//cout.setf(ios::scientific);

	if (globalArgs.fit == 0) {

		performRadius(globalArgs.output);

	} else {

		performFit(globalArgs.parameters,globalArgs.ffdata,globalArgs.csdata,globalArgs.output,globalArgs.fit);

	}
	
	/// End output

	time ( &rawtime );
	std::cout << "\n> End    : " <<  ctime (&rawtime) << "=== " << std::endl;
	
	return EXIT_SUCCESS;
}

// TODO:
