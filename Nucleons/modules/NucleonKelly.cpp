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

#include "NucleonKelly.h"

namespace ROOT {

	namespace Minuit2 {

// size = size of vector `a', e.~g., a0, a1 -> size=2
FFactorK::FFactorK (std::size_t size):
		a(size),
		b(size+2),
		c(2)
{
	numberOfParameters = size + size+2 + 2;
};

void FFactorK::LoadParameters (char* ds)
{
	ifstream myDataFile (ds);
	int por;
	int num = 0;
	double x;
	std::vector<double> fa, fb;
	if (myDataFile.is_open()) {
		while (myDataFile.peek() != EOF) {
			firstChar = myDataFile.peek();
			if ( (firstChar == '%') || (firstChar == '#') ) {
				getline (myDataFile,line);
			}
			else {
				myDataFile >> por >> x;
				if (num < a.size()) {
					a[por] = x;
				}
				if ((num >= a.size()) && (num < 2*a.size()+2)) {
					b[por] = x;
				}
				if ((num >= 2*a.size()+2)) {
					c[por] = x;
				}
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

void FFactorK::PrintParameters ()
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
	std::cout << "b[i]:" << std::endl;
	for (int i = 0; i < b.size(); ++i) {
		std::cout.width(6);
		std::cout << i;
		std::cout.width(10);
		std::cout << this->b[i];
		std::cout << std::endl;	
	}
	std::cout << "A, B:" << std::endl;
	for (int i = 0; i < c.size(); ++i) {
		std::cout.width(6);
		std::cout << i;
		std::cout.width(10);
		std::cout << this->c[i];
		std::cout << std::endl;	
	}
}

double FFactorK::G (double t, double M)
{
	double Qs = -t;
	double tau = Qs/(4.*M*M);
	double nom = 0.;
	double denom = 0.;

	for (int i = 0; i < a.size(); ++i) {
		nom = nom + a[i]*TMath::Power(tau,i);
		denom = denom + b[i]*TMath::Power(tau,i);
	}
	return nom/denom;
}

double FFactorK::GEP (double t)
{
	return FFactorK::G(t,massP);
}

double FFactorK::GMP (double t)
{
	return FFactorK::G(t,massP);
}

double FFactorK::GEN (double t, double M)
{
	double r;
	double Qs = -t;
	double tau = Qs/(4.*M*M);
	r = c[0]*tau/(1. + c[1]*tau)*this->GD(t);
	return r;
}

double FFactorK::GD (double t)
{
	double Qs = -t;
	double D = 1.+Qs/0.71;
	return 1./(D*D);
}

	}  // namespace Minuit2

}  // namespace ROOT
