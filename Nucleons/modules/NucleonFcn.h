/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonFcn.h $
 * $Id$
 *
 * @file
 * @brief	Header for FCN function of nucleon FF.
 */

#ifndef _NucleonFcn_H_
#define _NucleonFcn_H_

#include "Minuit2/FCNBase.h"


namespace ROOT {

	namespace Minuit2 {

class NucleonFcn : public FCNBase {

  private:
	const int fNPar;
	std::vector<int> fType;
	std::vector<double> fX;
	std::vector<double> fVal;
	std::vector<double> fErrUp;
	std::vector<double> fErrDown;
	std::vector<double> fTheta;
	std::vector<double> fEnergy;
	double fErrDef;

  public:
		
	NucleonFcn(
		const int nPar,
		const std::vector<int>& type,
		const std::vector<double>& x,
		const std::vector<double>& val,
		const std::vector<double>& errUp,
		const std::vector<double>& errDown,
		const std::vector<double>& theta,
		const std::vector<double>& energy):
		fNPar(nPar),
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

	~NucleonFcn () {}

};

	}  // namespace Minuit2
	
}  // namespace ROOT

#endif // _NucleonFcn_H_
