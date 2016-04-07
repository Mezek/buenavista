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

#include "MesonTest.h"

namespace ROOT {

	namespace Minuit2 {

FFactorE::FFactorE ()
{ 
	std::cout << ">> FFactor without defined size of parameters." << std::endl;

}

FFactorE::FFactorE ( std::size_t size ): a(size), v(size)
{ 
	t(1.,0.00001);
	std::cout << ">> FFactor is empty! Value of t: " << t << std::endl;

}

TComplex FFactorE::FFVal ( TComplex c )
{
	return c;
}

void FFactorE::testFunction ()
{
	std::cout << "Test passed" << std::endl;
}

	}  // namespace Minuit2

}  // namespace ROOT
