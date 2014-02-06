/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
/*
 * $Date: 2013-10-02 13:45:42 +0200 (Wed, 02 Oct 2013) $
 * $Revision: 383 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/Komparision.cpp $
 * $Id: Komparision.cpp 383 2013-10-02 11:45:42Z bartos $
 *
 * Description: Comparision of parametrizations
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
#include "TMath.h"
#include "TComplex.h"
#include "TStopwatch.h"

// 'A' = protonElectric
// 'B' = protonMagnetic
// 'C' = neutronElectric
// 'D' = neutronMagnetic

char dataFile[] = "dataNucleons.dat";
char dataFile0[] = "dataNucleonsCleared.dat";
//char dataFile1[] = "dataNucleonsApprox.dat";

char dataFileA[] = "../Data/dataProtonElectric.dat";
char dataFileB[] = "../Data/dataProtonMagnetic.dat";
char dataFileC[] = "../Data/dataNeutronElectric.dat";
char dataFileD[] = "../Data/dataNeutronMagnetic.dat";
char dataFileE[] = "../Data/dataProtonRatios.dat";
char dataFileF[] = "../Data/dataNeutronRatios.dat";

char parametersFile[] = "parNucleons-Kelly.dat";
char outputFile[] = "outNucleons-Kelly.dat";

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons.cpp"
#include "modules/ExperimentalData.cpp"
#include "modules/PlotGraph.cpp"
#include "modules/NucleonKelly.cpp"
//#include "NucleonFcn.cpp"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnPrint.h"

using namespace ROOT::Minuit2;

int main ( int argc, char **argv ) {

	TApplication theApp("Komparision", &argc, argv);

	time_t rawtime;
	time ( &rawtime );
	
	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2012 Erik Bartoš ===" << std::endl;

	// Start
	
	std::cout << "\n> Start  : " << ctime (&rawtime) << std::endl;
	
	// Set precision
	std::cout.precision(15); 

	std::cout << "\n> Fitted data:        '" << dataFile0 << "'" << std::endl;	
	std::cout << "> Start parameters:   '" << parametersFile << "'" << std::endl;
	std::cout << "> Output parameters:  '" << outputFile << "'" << std::endl;

	FFactorK pFF(4,4);
	pFF.LoadParameters(parametersFile);
	pFF.PrintParameters();
	
	int nPar = pFF.numberOfParameters;
	std::cout << "Number of FF parameters:       " << nPar << std::endl;
	//pFF.PrintParameters();

	// Get data
	
	ExperimentalData Z;
	Z.ReadData(dataFile0,3);
	int nData = Z.size();
	std::cout << "Number of fitted points:       " << nData << std::endl;
	Z.CheckData();

	std::vector<double> x = Z.X();
	std::vector<double> val = Z.Val();
	std::vector<double> errUp = Z.ErrUp();
	std::vector<double> errDown = Z.ErrDown();
	std::vector<int> type = Z.Type();
	
	// Create FCN function
	//NucleonFcn fFCN(nPar, Z.type, Z.x, Z.val, Z.errUp, Z.errDown);
	
	// Plot FFs with data points
	const int nPoints = 2500;
	double plotX[nPoints], plotY0[nPoints], plotY1[nPoints], plotY2[nPoints], plotY3[nPoints];
	double dataEX[nData], dataEY1[nData], dataEY2[nData];
	double tMin = -3.0;
	double tMax = 0.0;
	double tStep = (tMax-tMin)/nPoints;
	double tA = tMin;
	double res;
 	double gA,gB,gC,gD;
	
	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotX[i] = -tA;
		plotY0[i] = pFF.GEp(tA)/pFF.GD(tA);
		//std::cout <<  i << ": " << plotX[i] << "  " << plotY0[i] << std::endl;
    }

	// Experimental data
	double dip;
	for (int i = 0; i < nData; ++i) {
		dataEX[i] = -x[i];
		dip = 1.-x[i]/0.71;
		dataEY1[i] = val[i]*dip*dip;
		dataEY2[i] = val[i]/pFF.GD(x[i]);

		//std::cout <<  i << ": " << type[i] << "  " << x[i] << "  " << val[i] << std::endl;
		//std::cout <<  i << ": " << dataEX[i] << "  " << dataEY1[i] << "  " << dataEY2[i] << std::endl;
		//gC = 1./(dip*dip);
		//gD = pFF.GDipole(q2[i]);
		//std::cout <<  i << ": " << gC << "  " << gD << std::endl;
    }	

	PlotGraph graf1(10);
	const Char_t* title = "Kelly";
	graf1.view(nPoints,plotX,plotY0,title);

	//PlotGraph graf2;
	//graf2.viewPlusData(nPoints,plotX,plotY0,nData,dataEX,dataEY1,"Data");
	//graf2.viewData(nData,dataEX,dataEY1);
	
/*	// Create FCN function
	NucleonFcn fFCN(Z.q2, Z.val, Z.errUp, Z.errDown, Z.type);

	// Minuit
	double step = 0.01;
	MnUserParameters upar;
	upar.Add("t_1s  ", pFF.val[0], step, 2., 4.0);
	upar.Add("t_1v  ", pFF.val[1], step, 1.0, 4.0);
	upar.Add("t_2s  ", pFF.val[2], step, 1.7, 7.0);
	upar.Add("t_2v  ", pFF.val[3], step, 2.1, 7.5);
	upar.Add("f_Om  ", pFF.val[4], step, 0.1, 2.5);
	upar.Add("f_Ph  ", pFF.val[5], step, -0.2, 1.0);
	upar.Add("f_sOm ", pFF.val[6], step, -0.1, 1.0);
	upar.Add("f_sPh ", pFF.val[7], step, -0.8, 1.0);
	upar.Add("f_Rh  ", pFF.val[8], step, -0.1, 2.0);
	upar.Add("f_sRh ", pFF.val[9], step, -0.1, 3.0);
	upar.Add("f_OmM ", pFF.val[10], step, -1.0, 2.0);
	upar.Add("f_PhM ", pFF.val[11], step, -1.5, 4.7);
	upar.Add("f_sOmM", pFF.val[12], step, -2.5, 4.0);
	upar.Add("f_sPhM", pFF.val[13], step, -4.5, 4.0);
	upar.Add("f_RhM ", pFF.val[14], step, -1., 4.7);
	upar.Add("f_sRhM", pFF.val[15], step, -2.5, 4.5);

	// Create Migrad minimizer
	MnMigrad migrad(fFCN, upar, 1);

	/*migrad.Fix(0.);
	migrad.Fix(1.);
	migrad.Fix(2.);
	migrad.Fix(3.);
	migrad.Fix(4.);
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
	
	// Create Simplex minimizer
	//~MnSimplex simplex(fFCN, upar);

	// Minimize: migrad (maxfcn, tolerance)
/*
	timer.Start();
	//std::cout << "> ...1. minimalization..." << std::endl;	
	//FunctionMinimum min = migrad(100000,1.);

	//~MnMinos Minos(fFCN, min);

	//.migrad.Fix(15.);
	//std::cout << "...2. minimalization..." << std::endl;
	//FunctionMinimum min2 = migrad(10000,1.);
	
	std::cout << "> ...Last minimalization..." << std::endl;
	FunctionMinimum min9 = migrad(10000,1.);
	timer.Stop();
	
	// Standard output
	std::cout << "> Minimum: " << min9 << std::endl;
	std::cout << "> FCN value: " << min9.Fval() << std::endl;
	std::cout << "> Real time of minimalization: " << timer.RealTime() << " sec." << std::endl;

	// Scan and Plot parameters
 	/*{
    MnScan scan(fFCN, upar, 1);
    std::cout << "Scan parameters: " << scan.Parameters() << std::endl;
    MnPlot plot;
    for (unsigned int i = 0; i < upar.VariableParameters(); i++) {
		std::vector<std::pair<double, double> > xy = scan.Scan(i);
		//std::vector<std::pair<double, double> > xy = scan.scan(0);
		std::cout << i+1 << ". parameter - `" << upar.Name(i) << "'"<< std::endl;
		plot(xy);
	}
    std::cout << scan.Parameters() << std::endl;
	}*/

/*	// Output parameters
	ofstream myOutputParam (outputFile);
	for (int i = 0; i < nPar; i++) {
		myOutputParam << min9.UserParameters().Value(i) << std::endl;
		//std::cout << min9.UserParameters().Value(i) << std::endl;
	}
	myOutputParam.close();

	//nPlot.setParameters(min9.UserParameters());
*/
	
/*	// Plot FFs with data points
	const int nPoints = 2500;
	double plotX[nPoints], plotY0[nPoints], plotY1[nPoints], plotY2[nPoints], plotY3[nPoints];
	double dataX0[num0], dataX1[num1], dataX2[num2], dataX3[num3];
	double dataY0[num0], dataY1[num1], dataY2[num2], dataY3[num3];
	std::fill(dataX3, dataX3+num3, 0);
	double tMin = -4.5;
	double tMax = 15.0;
	double shag = (tMax-tMin)/nPoints;
	double tA = tMin;
	double res;
	TComplex z0(tMin,0.000001);
	TComplex z,zY,zW,gA,gB,gC,gD;

	
	for (int i = 0; i < nPoints; i++) {
		tA = tMin + i*shag;
		plotX[i] = tA;
 		//if (tA < t0v) { 
 		z = tMin + i*shag;
 		//}
 		//else {
 			//z = z0 + i*shag;
 		//}
		gA = nPlot.scalarE(z);
		gB = nPlot.vectorE(z);
		gC = nPlot.scalarM(z);
		gD = nPlot.vectorM(z);

		plotY0[i] = TComplex::Abs(gA + gB);
		plotY1[i] = TComplex::Abs(gC + gD);
		plotY2[i] = TComplex::Abs(gA - gB);
		plotY3[i] = TComplex::Abs(gC - gD);
    }
	
	int k = 0 , l = 0, m =0 , n = 0;
	for (int i = 0; i < q2.size(); i++) {
		switch ( type[i] ) {
			case 0:
				dataX0[k] = q2[i];
				dataY0[k] = val[i];
				//std::cout << ">> " << k << " " << type[i] << " " << dataX0[i] << " " << dataY0[k] <<std::endl;
				k++;
				break;
			case 1:
				dataX1[l] = q2[i];
				dataY1[l] = val[i];
				//std::cout << ">> " << l << " " << type[i] << " " << dataX1[i] << " " << dataY1[n] <<std::endl;
				l++;
				break;
			case 2:
				dataX2[m] = q2[i];
				dataY2[m] = val[i];
				//std::cout << ">> " << m << " " << type[i] << " " << dataX2[i] << " " << dataY2[n] <<std::endl;
				m++;
				break;
			case 3:
				dataX3[n] = q2[i];
				dataY3[n] = val[i];
				//std::cout << ">> " << n << " " << type[i] << " " << dataX3[i] << " " << dataY3[n] <<std::endl;
				n++;
				break;
			default:
				std::cout << "> Error: Unkonw type of data!" << std::endl;
		}
	}*/

/*	// Chi squared
	double thVal;
	double chi = 0.;
	double chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;
	TComplex cThVal;
	for (int i = 0; i < q2.size()-1; i++) {
		plotDx[i] = q2[i];
		plotDy[i] = val[i];
		//std::cout << i << " " << q2[i] << " " << val[i] << std::endl;
		//std::cout << i << " " << plotDx[i] << " " << plotDy[i] << std::endl;

		/*cThVal = nPlot.scalarM(q2[i])+nPlot.vectorM(q2[i]);
		thVal = cThVal.Re()*cThVal.Re()+cThVal.Im()*cThVal.Im();
		thVal = sqrt(thVal);

		delta = thVal-val[i];
		sigmaAverage = (errUp[i] + errDown[i])/2.;
		chi = delta/sigmaAverage;
		//std::cout << i << " " << thVal << " " << val[i] << " " << (nPlot.scalarM(q2[i])+nPlot.vectorM(q2[i]))*(nPlot.scalarM(q2[i])+nPlot.vectorM(q2[i])) << std::endl;
		chi2 += chi*chi;*/
		//plotDx[0] = q2[0];
		//std::cout << ">>" << q2.size() << std::endl;
	//}
	//std::cout << "Value of Ch2: " << chi2 << std::endl;
	//std::cout << "Value of Ch2/nData: " << chi2/nData << std::endl;*/

	/*for (int i = 0; i < q2.size(); i++) {
		std::cout << i << " " << plotDx[i] << " " << plotDy[i] << std::endl;
	}*/
	
	//std::cout << "t0s: " << t0s << std::endl;
	//std::cout << "t0v: " << t0v << std::endl;
	//std::cout << ">>> " << plotDx[0] << " " << plotDy[0] << std::endl;
	//std::cout << ">>> " << q2[0] << " " << val[0] << std::endl;

	// End output

	time ( &rawtime );
	std::cout << "\n> End    : " <<  ctime (&rawtime) << std::endl;
	theApp.Run(); //delete theApp;
	
	return EXIT_SUCCESS;
}

// To-do:
// - 
