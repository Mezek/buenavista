/**
 * $Date: 2013-10-08 11:35:13 +0200 (Tue, 08 Oct 2013) $
 * $Revision: 388 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonIachello.h $
 * $Id: NucleonIachello.h 388 2013-10-08 09:35:13Z bartos $
 *
 * @file
 * @brief	Header for Iachello form factor for nucleons.
 */

#ifndef NucleonUam_H_
#define NucleonUam_H_

#include "TMatrixDSym.h"
#include "TDecompChol.h"
#include "TRandom3.h"

namespace ROOT {

	namespace Minuit2 {

struct izo {
	izo() : name(), nor(), mesons(), tresh() {} // automatic initialization with constructor
	std::string name;
	TComplex nor;
	int mesons;
	int tresh;
};

struct hod {
	hod() : name(), val(), err(), down(), up() {} // automatic initialization with constructor
	//hod (std::string const& Name, double Val, double Down, double Up) : name(Name), val(Val), down(Down), up(Up) {}
	std::string name;
	double val, err, down, up;
};

class FFactor {

  private:
	double t;
	TComplex mS[6],wS[6],mV[6],wV[6],mwS2[6],mwV2[6];
	TComplex vM[6],vMc[6],mul[6],sub[6];
	izo FF[4];
	std::vector<hod> a;
	TMatrixD cov;
	int modelPar;
	bool handSome;
	
  public:
	FFactor ( std::size_t );
	int numberOfParameters;
	std::vector<hod> v;
    void LoadParameters ( char* );
	void SetParameter ( const int , const double );
    void SetParameters ( const std::vector<double>& );
	void FixParameters ( void );
    void PrintParameters ( void );
	void LoadCovMatrix ( char* );
	double A ( int );
	double E ( int );
	double Gfun (double t);
	double Afun (double t);
	double Ttrans (double t);
	double ScalarOne ( double );
	double ScalarTwo ( double );	
	double VectorOne ( double );
	double VectorTwo ( double );
	double GEP ( double );
	double GMP ( double );
	double GEN ( double );
	double GMN ( double );
	TComplex TypeDefVal ( const int &, const TComplex &, const double &, const double & );
	double SigmaTotalP ( const TComplex & );
	double SigmaTotalN ( const TComplex & );
	double Derive ( const TComplex &, const double );
	double DeriveXA ( const TComplex &, int, const double );
	double RadiusEP ( const double );
	double RadiusEPUncer ( const double );
	void RadiusEPUncerMC ( const double, const double );
	void CheckFormFactor ( const char*, double );
	void CheckParameters ( void );
	void CheckCovMatrix ( void );

	~FFactor () {};
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // NucleonUam_H_
