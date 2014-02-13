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

#ifndef _NucleonKelly_H_
#define _NucleonKelly_H_

namespace ROOT {

	namespace Minuit2 {
		
class FFactorK {

  private:
	std::vector<double> a, b;
	const int modelPar;
	double t;
	
  public:
	int numberOfParameters;
	//FFactorK () {}; // default constructor
	FFactorK ( std::size_t, int ); // for constructor initialization list
	void LoadParameters ( char* );
	void PrintParameters ();
	double G ( double, double );
	double GEp ( double );
	double GMp ( double );
	double GEn ( double );
	double GMn ( double );
	double GD ( double );
	~FFactorK () {};
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // _NucleonKelly_H_
