/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Header for FCN function of transition FF.
 */

#ifndef _TransitionFcn_H_
#define _TransitionFcn_H_

#include "Minuit2/FCNBase.h"


namespace ROOT {

	namespace Minuit2 {

class TransitionFcn : public FCNBase {

  private:
	const int fTParticle;
	std::vector<int> fType;
	std::vector<double> fX;
	std::vector<double> fVal;
	std::vector<double> fErrUp;
	std::vector<double> fErrDown;
	std::vector<double> fTheta;
	std::vector<double> fEnergy;
	double fErrDef;

  public:
		
	TransitionFcn(
		const int tParticle,
		const std::vector<int>& type,
		const std::vector<double>& x,
		const std::vector<double>& val,
		const std::vector<double>& errUp,
		const std::vector<double>& errDown,
		const std::vector<double>& theta,
		const std::vector<double>& energy):
		fTParticle(tParticle),
		fType(type),
		fX(x),
		fVal(val),
		fErrUp(errUp),
		fErrDown(errDown),
		fTheta(theta),
		fEnergy(energy),
		fErrDef(1.) {}

	virtual double Up() const {return fErrDef;}
	virtual double operator() ( const std::vector<double>& ) const;

	std::vector<int> Type() const { return fType; }
	std::vector<double> X() const { return fX; }
	std::vector<double> Val() const { return fVal; }
	std::vector<double> ErrUp() const { return fErrUp; }
	std::vector<double> ErrDown() const { return fErrDown; }
	std::vector<double> Theta() const { return fTheta; }
	std::vector<double> Energy() const { return fEnergy; }

	void SetErrorDef(double def) { fErrDef = def; }

	~TransitionFcn () {}

};

	}  // namespace Minuit2
	
}  // namespace ROOT

#endif // _TransitionFcn_H_
