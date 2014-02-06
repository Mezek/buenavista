/**
 * $Date: 2013-05-13 15:39:22 +0200 (Mon, 13 May 2013) $
 * $Revision: 331 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonKelly.cpp $
 * $Id: NucleonKelly.cpp 331 2013-05-13 13:39:22Z bartos $
 *
 * @file
 * @brief	FCN function of Kelly's parametrization.
 */

#include "NucleonKelly.h"

namespace ROOT {

	namespace Minuit2 {

FFactorK::FFactorK (std::size_t size, int poke):
		a(size),
		b(size),
		modelPar(poke)
{
	numberOfParameters = 2*modelPar;
};

void FFactorK::LoadParameters (char* ds)
{
	ifstream myDataFile (ds);
	double x,y;
	int num = 0;
	if (myDataFile.is_open()) {
		while (myDataFile.peek() != EOF) {
			myDataFile >> a[num] >> b[num];
			++num;
		}
		myDataFile.close();
		num = num-1;
		if (num != a.size()) {
			std::cout << ">> Error: Number of parameters in " << ds << ": " << num << "' differs from declared: " << a.size() << "!" << std::endl;
		} else {
			std::cout << "\n " << std::endl;
		}
		/*for (int i = 0; i < a.size(); ++i) {
			std::cout << i << ": " << a[i] << " " << b[i] << std::endl;
		}*/
	}
	else std::cerr << ">> Error: Unable to open parametric file: '" << ds << "'!" << std::endl;
}

void FFactorK::PrintParameters ()
{
	std::cout << "\n>> Actual pararameters:" << std::endl;
	std::cout << "\ti" << "  \t" << "a[i]" << "\t " << "b[i]" << std::endl;
	std::cout << "    -----------------------------" << std::endl;
	for (int i=0; i<a.size(); ++i) {
		//std::cout << i+1 << ". " << a[i] << std::endl;
		std::cout << "\t" << i << ". \t" << this->a[i] << "\t " << this->b[i] << std::endl;
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

double FFactorK::GEp (double t)
{
	return FFactorK::G(t,massP);
}

double FFactorK::GMp (double t)
{
	return FFactorK::G(t,massP);
}

double FFactorK::GD (double t)
{
	double Qs = -t;
	double D = 1.+Qs/0.71;
	return 1./(D*D);
}

	}  // namespace Minuit2

}  // namespace ROOT
