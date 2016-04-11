/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonFcn.cpp $
 * $Id$
 *
 * @file
 * @brief	FCN function of transition FF.
 */

#include "TransitionFcn.h"

#include <cassert>

namespace ROOT {

	namespace Minuit2 {

double TransitionFcn::operator() (const std::vector<double>& par) const {
  
	TComplex nFv,nE,nM;
	double rate;

	FFactorT bornTBW(fTParticle);

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
		nFv = bornTBW.FFAbsVal(fX[n]);
		delta = nFv.Re()-fVal[n];
		sigmaAverage = (fErrUp[n] + fErrDown[n])/2.;
		chi = delta/sigmaAverage;
		chi2 = chi*chi;
		Chi2 += chi2;
	}

	return Chi2;
}

	}  // namespace Minuit2

}  // namespace ROOT
