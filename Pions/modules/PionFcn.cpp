/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	FCN function of pion FF.
 */

#include "PionFcn.h"

#include <cassert>

namespace ROOT {

	namespace Minuit2 {

double NucleonFcn::operator() (const std::vector<double>& par) const {
  
	//assert(par.size() == modelPar);

	FFactor pionik(fNPar);

	pionik.SetParameters(par);
	//pionik.PrintParameters();
	
	double chi = 0.;
	double chi2 = 0.;
	double Chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;

	TComplex PV;
	TComplex pFv;

	for(unsigned int n = 0; n < fX.size(); ++n) {
		/**  default type */
		if ((fType[n] == 0) || (fType[n] == 1)) {
			PV = pionik.Value(fX[n]);
			pFv = pFv.Abs(PV);
		}
		if (fType[n] >= 2) {
			std::cout << ">> Warning! No type of form factor defined!" << std::endl;
		}

		delta = pFv.Re()-fVal[n];
		sigmaAverage = (fErrUp[n] + fErrDown[n])/2.;
		chi = delta/sigmaAverage;
		chi2 = chi*chi;

		Chi2 += chi2;
	}

	return Chi2;
}

	}  // namespace Minuit2

}  // namespace ROOT
