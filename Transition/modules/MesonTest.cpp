/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonUam.cpp $
 * $Id$
 *
 * @file
 * @brief	Model of U&A transition form factors: Eta, EtaPrime, PiZero.
 */

#include "MesonUam.h"

namespace ROOT {

	namespace Minuit2 {
		
FFactor::FFactor ()
{ 
	t(1.,0.00001);
}

void FFactor::testFunction (double a)
{
	std::cout << "Test a=" << a << std::endl;
	std::cout << "Test t=" << t << std::endl;
}

	}  // namespace Minuit2

}  // namespace ROOT
