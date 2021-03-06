/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonUam.h $
 * $Id$
 *
 * @file
 * @brief	Header for U&A form factor for nucleons.
 */

#ifndef _NucleonUam_H_
#define _NucleonUam_H_

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
	TComplex t;
	TComplex mS[6],wS[6],mV[6],wV[6],mwS2[6],mwV2[6];
	TComplex vM[6],vMc[6],mul[6],sub[6];
	izo FF[4];
	std::vector<hod> a;
	std::string FFtype[4];
	TMatrixD cov;
	int modelPar;
	int expressPar;
	int allMesons;
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
	void ExpressedParameters ( void );
	void LoadCovMatrix ( char* );
	double A ( int );
	double E ( int );
	TComplex W ( const TComplex &, const TComplex &, const TComplex &, const double );
	TComplex ScalarOne ( TComplex );
	TComplex ScalarTwo ( TComplex );	
	TComplex VectorOne ( TComplex );
	TComplex VectorTwo ( TComplex );
	TComplex GEP ( TComplex );
	TComplex GMP ( TComplex );
	TComplex GEN ( TComplex );
	TComplex GMN ( TComplex );
	TComplex TypeDefVal ( const int &, const TComplex &, const double &, const double & );
	double AbsGEP ( TComplex );
	double AbsGMP ( TComplex );
	double AbsGEN ( TComplex );
	double AbsGMN ( TComplex );
	double SigmaTotalP ( const TComplex & );
	double SigmaTotalN ( const TComplex & );
	double Derive ( const int, const TComplex &, const double );
	double DeriveXA ( const int, const TComplex &, int, const double );
	double Radius ( const int, const double );
	double RadiusUncer ( const int, const double );
	void RadiusUncerMC ( const int, const double, const double );
	double DeriveOld ( const TComplex &, const double );
	double DeriveXAOld ( const TComplex &, int, const double );
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

#endif // _NucleonUam_H_
