/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Header for arbitrary UA form factor of pion + Omega-Rho interference.
 */

#ifndef _PionUam_H_
#define _PionUam_H_

namespace ROOT {

	namespace Minuit2 {

struct hod {
	hod() : name(), val(), err(), down(), up() {} // automatic initialization with constructor
	std::string name;
	double val, err, down, up;
};

class FFactor {

  private:
	TComplex t;
	std::vector<hod> a;
	TMatrixD cov;
	int modelPar;
	bool handsome;

  public:
	FFactor ( std::size_t );
	int numberOfParameters;
	void LoadParameters ( char* );
	void SetParameter ( const int, const double );
	void SetParameters ( const std::vector<double>& );
	void PrintParameters ( void );
	double A ( int );
	double E ( int );
	void LoadCovMatrix ( char* );
	TComplex W ( const TComplex &, const double);
	TComplex Value ( TComplex );
	double ValueSquared ( TComplex );

	~FFactor () {};
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // _PionUam_H_
