/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/modules/NucleonUamH.h $
 * $Id$
 *
 * @file
 * @brief	Header for U&A form factor for nucleons, v. hyperon.
 */

#ifndef NucleonUamH_H_
#define NucleonUamH_H_

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
	TMatrixD cov;
	int modelPar;
	bool handSome;
	int particleType;
	double mS2[6],mV2[6];
	
  public:
	double massH1, massH2;
	std::string particleName1, particleName2;
	FFactor ();
	FFactor ( int );
	int numberOfParameters;
	std::vector<hod> v;
	void SetParticle ( int );
	int GetParticle ( void );
	void LoadParameters ( char* );
	void SetParameter ( const int , const double );
	void SetParameters ( const std::vector<double>& );
	void FixParameters ( void );
	void PrintParameters ( void );
	void LoadCovMatrix ( char* );
	double A ( int );
	double E ( int );
	void TransformSU3 ( void );
	TComplex W ( const TComplex &, const TComplex &, const TComplex &, const double );
	TComplex ScalarOne ( TComplex );
	TComplex ScalarTwo ( TComplex );	
	TComplex VectorOne ( TComplex );
	TComplex VectorTwo ( TComplex );
	TComplex GE1 ( TComplex );
	TComplex GM1 ( TComplex );
	TComplex GE2 ( TComplex );
	TComplex GM2 ( TComplex );
	TComplex TypeDefVal ( const int &, const TComplex &, const double &, const double & );
	double GE1N ( void );
	double GM1N ( void );
	double GE2N ( void );
	double GM2N ( void );
	double AbsGE1 ( TComplex );
	double AbsGM1 ( TComplex );
	double AbsGE2 ( TComplex );
	double AbsGM2 ( TComplex );
	double SigmaTotal1 ( const TComplex & );
	double SigmaTotal2 ( const TComplex & );
	double RadiusE1 ( const double );
	double RadiusE2 ( const double );
	double RadiusM1 ( const double );
	double RadiusM2 ( const double );
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

#endif // NucleonUamH_H_
