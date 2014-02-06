/**
 * $Date: 2013-05-13 15:39:22 +0200 (Mon, 13 May 2013) $
 * $Revision: 331 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonKelly.h $
 * $Id: NucleonKelly.h 331 2013-05-13 13:39:22Z bartos $
 *
 * @file
 * @brief	Header for FCN function of Kelly's parametrization.
 */

#ifndef NucleonKelly_H_
#define NucleonKelly_H_

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

#endif // NucleonKelly_H_
