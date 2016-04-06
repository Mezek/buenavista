/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonUam.h $
 * $Id$
 *
 * @file
 * @brief	Header for U&A transition form factors: Eta, EtaPrime, PiZero.
 */

#ifndef _MesonTest_H_
#define _MesonTest_H_

#include "TMatrixDSym.h"
#include "TDecompChol.h"
#include "TRandom3.h"

namespace ROOT {

	namespace Minuit2 {

class FFactor {

  private:
	TComplex t;
	std::vector<hod> a;
	std::string FFtype[4];
	
  public:
	FFactor ();
	double testFunction ( const double );

	~FFactor () {};
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // _MesonTest_H_
