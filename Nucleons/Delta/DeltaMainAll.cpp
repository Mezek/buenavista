/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Main programe for delta nucleon form factors.
 *  
 * <b>Compilation:</b>
 * @code
 * > g++ -o DeltaMain.exe DeltaMain.cpp `root-config --cflags` `root-config --libs` -lMinuit2
 * @endcode
 * 
 * <b>Usage</b> (to run as standalone application in console):
 * @code
 * > ./DeltaMain.exe --help
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
/** Transition data files.*/
char dataFileA[] = "../../Data/dataGMstar.dat";
char dataFileB[] = "../../Data/dataREM.dat";
char dataFileC[] = "../../Data/dataRSM.dat";
///@}

char parametersFile[] = "parElasticFF.dat";            ///< Input parameters.
char outputFile[] = "outDelta.dat";                    ///< Output file.

#include "../modules/ConstBasic.cpp"
#include "../modules/ConstMesons.cpp"
#include "../modules/ConstDelta.cpp"
#include "../modules/Nucleon3G-D.cpp"
#include "../modules/NucleonFcnCS.cpp"

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

#include "DeltaData.cpp"
#include "DeltaPlot.cpp"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	time_t rawtime;
	time ( &rawtime );
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2015 Erik Bartoš" << std::endl;

	/// Parse options

	//performUsage(argc, argv);

	/// Start
	
	TApplication theApp("DeltaNucleon", &argc, argv);
	std::cout << "\n> Start  : " << ctime (&rawtime) << std::endl;
	
	/// Set precision
	std::cout.precision(15); 
	//cout.setf(ios::scientific);

	int graph = 2;
	if (argc > 1 ) { graph = atoi(argv[1]); }
	std::string dataFile;

	switch (graph) {
		case 0: {

			/// Plot G_M^*
			dataFile = "../../Data/dataGMstar.dat";
				
			break;
		}
		case 1: {

			/// Plot R_EM
			dataFile = "../../Data/dataREM.dat";
				
			break;
		}
		case 2: {

			/// Plot R_SM
			dataFile = "../../Data/dataRSM.dat";
				
			break;
		}
		default: 
			break;
	}
	DeltaData A;
	A.ReadData(dataFile.c_str());	
	A.ShowData();

	int num = A.size();
	std::cout << "\nNumber of plotted points:      " << num << std::endl;
	std::vector<double> x = A.X();
	std::vector<double> val = A.Val();
	std::vector<double> errUp = A.ErrUp();
	std::vector<double> errDown = A.ErrDown();

	double X[num], Y[num], U[num], D[num];
	int k = 0;
	for (int i = 0; i < num; i++) {
		X[k] = -x[i];
		//Y[k] = val[i]/3.*(1.+X[k]/0.71)*(1.+X[k]/0.71);
		Y[k] = val[i];
		D[k] = errDown[i];
		U[k] = errUp[i];
		++k;
	}	
	DeltaPlot graf(2);
	graf.viewData(num, X, Y, "Data");

	FFactor GMS(12);
	GMS.LoadParameters(parametersFile);
	GMS.CheckParameters();
	//GMS.PrintParameters();

	const int nPoints = 2500;
	double qMin;
	double qMax;
	double qStep;
	double qA;

	double plotGX[nPoints], plotGY[nPoints];
	qMin = 0.004;
	qMax = 3.0;
	qStep = (qMax-qMin)/nPoints;
	qA = qMin;
	
	for (int i = 0; i < nPoints; i++) {
		qA = qMin + i*qStep;
		double qA2 = qA*qA;
		double gen = GMS.AbsGEN(-qA2);
		double gmn = GMS.AbsGMN(-qA2);
		double msq = massDi*massDi - massNucl*massNucl - qA2;
		double abq = sqrt((qA2 + msq*msq)/(4.*massDi*massDi));
		plotGX[i] = qA2;
		//plotGY[i] = abq*massNucl/qA2*gen/gmn;
		//plotGY[i] = gen/gmn;
		//std::cout << massNucl << " " << gen/qA2 << std::endl;
		double mPlu = 0.5*(massDi + massNucl);
		double mNom = 1. - qA2/(4.*mPlu*mPlu);
		double mDen = 1. - qA2/(4.*massNucl*massNucl);
		//plotGY[i] = sqrt(2.)/3.*(massDi+massNucl)/massNucl*mNom/mDen*gmn;
		double mDip = 1. + qA2/0.71;
		double gD = 1./mDip/mDip;
		double gDmn = gD*(-muN);
		plotGY[i] = sqrt(2.)*2./3.*gDmn;
		//plotGY[i] = sqrt(2.)*2./3.*gmn;
	}

	graf.viewPlusData(nPoints, plotGX, plotGY, num, X, Y, "title");
	
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
