/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonFcn.cpp $
 * $Id$
 *
 * @file
 * @brief	FCN function of nucleon FF.
 */

#include "NucleonFcn.h"

#include <cassert>

namespace ROOT {

	namespace Minuit2 {

double NucleonFcn::operator() (const std::vector<double>& par) const {
  
	//assert(par.size() == modelPar);

	TComplex nFv,nE,nM;
	double rate;

	FFactor bornTBW(fNPar);

	bornTBW.SetParameters(par);
	//bornTBW.PrintParameters();
	
	double chi = 0.;
	double chi2 = 0.;
	double Chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;
	double Asym = 0.;
	TComplex Fs1,Fs2,Fv1,Fv2;

	for(unsigned int n = 0; n < fX.size(); ++n) {
		/**  "protonElectric" */
		if ((fType[n] == 1) || (fType[n] == 2)) {
			nFv = bornTBW.AbsGEP(fX[n]);
			//std::cout << n << " " << fType[n] << std::endl;
		}
		/**  "protonMagnetic" */
		if ((fType[n] == 3) || (fType[n] == 4)) {
			nFv = bornTBW.AbsGMP(fX[n]);
		}
		/**  "neutronElectric" */
		if ((fType[n] == 5) || (fType[n] == 6)) {
			nFv = bornTBW.AbsGEN(fX[n]);
		}
		/**  "neutronMagnetic" */
		if (fType[n] == 7) { nFv = bornTBW.GMN(fX[n]); }
		if (fType[n] == 8) { nFv = bornTBW.AbsGMN(fX[n]); }
		/**  "protonRatios" */
		if (fType[n] == 9) {
			nE = bornTBW.AbsGEP(fX[n]);
			nM = bornTBW.AbsGMP(fX[n]);
			nFv = nFv.Abs((1.+ammP)*nE/nM);
		}
		/**  "neutronRatios" */
		if (fType[n] == 10) {
			nE = bornTBW.AbsGEN(fX[n]);
			nM = bornTBW.AbsGMN(fX[n]);
			nFv = nFv.Abs(ammN*nE/nM);
		}
		if (fType[n] >= 11) {
			std::cout << ">> Warning! No type of form factor defined!" << std::endl;
		}

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
			chi2 = chi*chi-2.*Asym*chi*chi*chi+5.*Asym*Asym*chi*chi*chi*chi;
			//chi2 = chi*chi;
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
