/**
 * $Date: 2013-05-22 10:03:30 +0200 (Wed, 22 May 2013) $
 * $Revision: 348 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/Mesons.cpp $
 * $Id: Mesons.cpp 348 2013-05-22 08:03:30Z bartos $
 *
 * @file
 * @brief	Calculate some values for mesons.
 *  
 * <b>Compilation:</b>
 * @code
 * > g++ -o Mesons.exe Mesons.cpp `root-config --cflags` `root-config --libs`
 * @endcode
 */

#include <iostream>
#include <fstream>
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

#include "modules/ConstBasic.cpp"
#include "modules/ConstMesons.cpp"

using namespace std;

struct izo {
	izo() : name(), mass(), width(), widthEE() {} // automatic initialization with constructor
	std::string name;
	double mass;
	double width;
	double widthEE;
};

double uniF ( double m, double W ) {
	double s = 0.;
	s = 2.*alpha*sqrt(TMath::Pi()*m/3./W);
	return s;
}

int main ( int argc, char **argv ) {

	std::cout << "=== Program: " << argv[0] << std::endl;
	std::cout << "=== Version 1.0, (c) 2012 Erik Barto¨" << std::endl;

	izo scalar[6], vector[6];
	// Names
	scalar[0].name = "Om";
	scalar[1].name = "Phi";
	scalar[2].name = "Om1P";
	scalar[3].name = "Phi1P";
	scalar[4].name = "Om2P";
	scalar[5].name = "Phi2P";

	vector[0].name = "Rho";
	vector[1].name = "Rho1P";
	vector[2].name = "Rho2P";
	vector[3].name = "Rho3P";
	vector[4].name = "Rho4P";
	vector[5].name = "Rho5P";
		
	// Masses and widths
	scalar[0].mass = massOm;
	scalar[1].mass = massPhi;
	scalar[2].mass = massOm1P;
	scalar[3].mass = massPhi1P;
	scalar[4].mass = massOm2P;
	scalar[5].mass = massPhi2P;

	scalar[0].width = widthOm;
	scalar[0].widthEE = widthEEOm;
	scalar[1].width = widthPhi;
	scalar[1].widthEE = widthEEPhi;
	scalar[2].width = widthOm1P;
	scalar[2].widthEE = widthEEOm1P;
	scalar[3].width = widthPhi1P;
	scalar[4].width = widthOm2P;
	scalar[4].widthEE = widthEEOm2P;
	scalar[5].width = widthPhi2P;
	
	vector[0].mass = massRho;
	vector[1].mass = massRho1P;
	vector[2].mass = massRho2P;
	vector[3].mass = massRho3P;
	vector[4].mass = 0.;
	vector[5].mass = 0.;

	vector[0].width = widthRho;
	vector[0].widthEE = widthEERho;
	vector[1].width = widthRho1P;
	vector[1].widthEE = widthEERho1P;
	vector[2].width = widthRho2P;
	vector[2].widthEE = widthEERho2P;
	vector[3].width = widthRho3P;
	vector[4].width = 0.001;
	vector[5].width = 0.001;
	
	for (int i = 0; i < 6; ++i) {
		std::cout << scalar[i].name << "\t" << scalar[i].mass << "\t" << scalar[i].width << "\t" << std::endl;
	}
	std::cout << "\n";
	for (int i = 0; i < 6; ++i) {
		std::cout << vector[i].name << "\t" << vector[i].mass << "\t" << vector[i].width << "\t" << std::endl;
	}

	std::cout << "\n";
	double fz[10], fp[10], fd[10];
	fz[0] = uniF(scalar[0].mass,scalar[0].widthEE);
	fz[1] = uniF(scalar[1].mass,scalar[1].widthEE);
	fz[2] = uniF(vector[0].mass,vector[0].widthEE);
	
	std::cout << "f_" << scalar[0].name << "\t" << fz[0] << "\t" << fz[0]*fz[0]/4./TMath::Pi() << "\t" << std::endl;
	std::cout << "f_" << scalar[1].name << "\t" << fz[1] << "\t" << fz[1]*fz[1]/4./TMath::Pi() << "\t" << std::endl;
	std::cout << "f_" << vector[0].name << "\t" << fz[2] << "\t" << fz[2]*fz[2]/4./TMath::Pi() << "\t" << std::endl;

	double p[3];
	p[0] = fz[1]/fz[0];
	p[1] = fz[2]/fz[0];
	p[2] = fz[2]/fz[1];
	std::cout << p[0] << " : " << p[1] << " : " << p[2] << std::endl;
	
	fp[0] = uniF(scalar[2].mass,scalar[2].widthEE);
	fp[2] = uniF(vector[1].mass,vector[1].widthEE);
	std::cout << "\n";
	std::cout << "f_" << scalar[2].name << "\t" << fp[0] << "\t" << std::endl;
	std::cout << "f_" << vector[1].name << "\t" << fp[2] << "\t" << std::endl;
	fp[1] = fp[0]*p[0];
	std::cout << "C: f_" << scalar[3].name << "\t" << fp[1] << "\t" << std::endl;
	fp[1] = fp[2]/p[2];
	std::cout << "C: f_" << scalar[3].name << "\t" << fp[1] << "\t" << std::endl;

	fd[0] = uniF(scalar[4].mass,scalar[4].widthEE);
	fd[2] = uniF(vector[2].mass,vector[2].widthEE);
	std::cout << "\n";
	std::cout << "f_" << scalar[4].name << "\t" << fd[0] << "\t" << std::endl;
	std::cout << "f_" << vector[2].name << "\t" << fd[2] << "\t" << std::endl;
	fd[1] = fd[0]*p[0];
	std::cout << "C: f_" << scalar[5].name << "\t" << fd[1] << "\t" << std::endl;
	fd[1] = fd[2]/p[2];
	std::cout << "C: f_" << scalar[5].name << "\t" << fd[1] << "\t" << std::endl;

	//std::cout << u[0] << std::endl;
	
}
