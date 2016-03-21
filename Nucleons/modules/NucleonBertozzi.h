/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonKelly.h $
 * $Id$
 *
 * @file
 * @brief	Header for FCN function of Kelly's parametrization.
 */

#ifndef _NucleonBertozzi_H_
#define _NucleonBertozzi_H_

namespace ROOT {

	namespace Minuit2 {
		
class FFactorB {

  private:
	std::vector<double> a, b, c;
	double t;
	char firstChar;
	std::string line;
	
  public:
	int numberOfParameters;
	//FFactorK () {}; // default constructor
	FFactorB ( std::size_t ); // for constructor initialization list
	void LoadParameters ( char* );
	void PrintParameters ();
	double G ( double, double );
	double GEP ( double );
	double GMP ( double );
	double GEN ( double );
	double GMN ( double );
	double GD ( double );
	~FFactorB () {};
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // _NucleonBertozzi_H_
