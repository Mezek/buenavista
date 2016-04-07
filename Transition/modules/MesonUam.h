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

struct izo {
	izo() : name(), nor(), mesons() {} // automatic initialization with constructor
	std::string name;
	TComplex nor;
	int mesons;
};

struct hod {
	hod() : name(), val(), err(), down(), up() {} // automatic initialization with constructor
	std::string name;
	double val, err, down, up;
};

class FFactorT {

  private:
	TComplex t;
	TComplex mS[6],wS[6],mV[6],wV[6],mwS2[6],mwV2[6];
	TComplex vM[6],vMc[6],mul[6],sub[6];
	izo FF;
	std::vector<hod> a;
	TMatrixD cov;
	int modelPar;
	int allMesons;
	bool handSome;
	
  public:
	FFactorT ();
	FFactorT ( int );
	int numberOfParameters;
	std::vector<hod> v;
	void LoadParameters ( char* );
	void SetParameter ( const int , const double );
	void SetParameters ( const std::vector<double>& );
	void FixParameters ( void );
	void PrintParameters ( void );
	double A ( int );
	double E ( int );
	TComplex W ( const TComplex &, const TComplex &, const TComplex &, const double );
	TComplex ScalarP ( TComplex );
	TComplex VectorP ( TComplex );
	TComplex FFVal ( TComplex );
	double FFAbsVal ( TComplex );

	~FFactorT () {};
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // _MesonTest_H_
