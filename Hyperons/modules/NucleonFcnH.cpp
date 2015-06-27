/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/modules/NucleonFcnH.cpp $
 * $Id$
 *
 * @file
 * @brief	FCN function for nucleon FF, v. hyperon.
 */

#include "NucleonFcnH.h"

#include <cassert>

namespace ROOT {

	namespace Minuit2 {

double NucleonFcn::operator() (const std::vector<double>& par) const {
  
	//assert(par.size() == modelPar);

	FFactor bornTBW;

	bornTBW.SetParameters(par);
	//bornTBW.PrintParameters();
	
	double chi = 0.;
	double chi2 = 0.;
	double Chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;
	double Asym = 0.;
	TComplex nFv,nE,nM;
	TComplex Fs1,Fs2,Fv1,Fv2;

	for(unsigned int n = 0; n < fX.size(); ++n) {

		nFv = bornTBW.TypeDefVal(fType[n],fX[n],fTheta[n],fEnergy[n]);
		delta = nFv.Re()-fVal[n];
		sigmaAverage = (fErrUp[n] + fErrDown[n])/2.;
		chi = delta/sigmaAverage;
		chi2 = chi*chi;

		/*
		delta = nFv.Re()-fVal[n];
		sigmaAverage = (fErrUp[n] + fErrDown[n])/2.;
		Asym = (fErrUp[n] - fErrDown[n])/(fErrUp[n] + fErrDown[n]);
		chi = delta/sigmaAverage;
		if (Asym == 0.) { chi2 = chi*chi; } else {
			//chi2 = chi*chi-2.*Asym*chi*chi*chi+5.*Asym*Asym*chi*chi*chi*chi;
			chi2 = chi*chi;
		}
		*/

		Chi2 += chi2;

		//std::cout << ">> Function: " << n << " " << sigmaAverage << " " << chi << " " << chi2 << std::endl;
		//std::cout << ">> Function: " << n << " " << rate << " " << nFv << std::endl;
		//std::cout << ">> Function: " << n << " " << chi << " " << sigmaAverage << std::endl;
		//std::cout << ">> Function: " << n << " " << chi << " " << chi2 << std::endl;

		//std::cout << n << " " << fType[n] << " " << nFv << " " << fVal[n] << std::endl;
		//std::cout << ">> " << fX.size() << std::endl;
	}

	return Chi2;
}

	}  // namespace Minuit2

}  // namespace ROOT
