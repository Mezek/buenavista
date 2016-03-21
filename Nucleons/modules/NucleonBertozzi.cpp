/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonKelly.cpp $
 * $Id$
 *
 * @file
 * @brief	FCN function of Kelly's parametrization.
 */

#include "NucleonBertozzi.h"

namespace ROOT {

	namespace Minuit2 {

// size = size of vector `a', e.~g., a0, a1 -> size=2
FFactorB::FFactorB (std::size_t size):
		a(size)
{
	numberOfParameters = size;
};

void FFactorB::LoadParameters (char* ds)
{
	ifstream myDataFile (ds);
	int por;
	int num = 0;
	double x, errx;
	std::vector<double> fa, fb;
	if (myDataFile.is_open()) {
		while (myDataFile.peek() != EOF) {
			firstChar = myDataFile.peek();
			if ( (firstChar == '%') || (firstChar == '#') ) {
				getline (myDataFile,line);
			}
			else {
				myDataFile >> por >> x >> errx;
				a[por] = x;
				getline (myDataFile,line);
				++num;
			}
		}
		myDataFile.close();
		if (num != numberOfParameters) {
			std::cout << ">> Error: Number of parameters in `" << ds << "': "
				      << num << " differs from declared: " << numberOfParameters << "!" << std::endl;
		} else {
			//std::cout << "\n " << std::endl;
		}
	}
	else std::cerr << ">> Error: Unable to open parametric file: '" << ds << "'!" << std::endl;
}

void FFactorB::PrintParameters ()
{
	std::cout << "\n> Kelly's form factor:" << std::endl;
	std::cout << "> Actual pararameters:" << std::endl;
	std::cout << "a[i]:" << std::endl;
	for (int i = 0; i < a.size(); ++i) {
		std::cout.width(6);
		std::cout << i;
		std::cout.width(10);
		std::cout << this->a[i];
		std::cout << std::endl;	
	}
}

double FFactorB::G (double t, double M)
{
	double Qs = -t;
	double tau = Qs/(4.*M*M);
	double nom = 0.;
	double denom = 0.;

	for (int i = 0; i < a.size(); ++i) {
		nom += 1.;
		denom += 1.;
	}
	return nom/denom;
}

double FFactorB::GEP (double t)
{
	return FFactorB::G(t,massP);
}

double FFactorB::GMP (double t)
{
	return FFactorB::G(t,massP);
}

double FFactorB::GEN (double t)
{
	double r;
	double Qs = -t*1000000./hTransC/hTransC;
	double r12 = a[2]*a[2] + a[1]/(2.*a[0]);
	double r22 = a[2]*a[2] - a[1]/(2.*a[0]);
	r = a[0]/TMath::Power(1. + Qs*r12/12., 2.) - a[0]/TMath::Power(1. + Qs*r22/12., 2.);
	return r;
}

double FFactorB::GD (double t)
{
	double Qs = -t;
	double D = 1.+Qs/0.71;
	return 1./(D*D);
}

	}  // namespace Minuit2

}  // namespace ROOT
