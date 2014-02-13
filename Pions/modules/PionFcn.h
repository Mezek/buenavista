/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Header for FCN function of pion FF.
 */

#ifndef _PionFcn_H_
#define _PionFcn_H_

#include "Minuit2/FCNBase.h"


namespace ROOT {

	namespace Minuit2 {

class PionFcn : public FCNBase {

  private:
	const int fNPar;
	std::vector<int> fType;
	std::vector<double> fX;
	std::vector<double> fVal;
	std::vector<double> fErrUp;
	std::vector<double> fErrDown;
	double fErrDef;

  public:
		
	PionFcn(
		const int nPar,
		const std::vector<int>& type,
		const std::vector<double>& x,
		const std::vector<double>& val,
		const std::vector<double>& errUp,
		const std::vector<double>& errDown):
		fNPar(nPar),
		fType(type),
		fX(x),
		fVal(val),
		fErrUp(errUp),
		fErrDown(errDown),
		fErrDef(1.) {}

	virtual double Up() const {return fErrDef;}
	virtual double operator() ( const std::vector<double>& ) const;

	std::vector<int> Type() const { return fType; }
	std::vector<double> X() const { return fX; }
	std::vector<double> Val() const { return fVal; }
	std::vector<double> ErrUp() const { return fErrUp; }
	std::vector<double> ErrDown() const { return fErrDown; }

	void SetErrorDef(double def) { fErrDef = def; }

	~PionFcn () {}

};

	}  // namespace Minuit2
	
}  // namespace ROOT

#endif //_PionFcn_H_
