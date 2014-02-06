/**
 * $Date: 2013-05-13 15:39:22 +0200 (Mon, 13 May 2013) $
 * $Revision: 331 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonFcn.h $
 * $Id: NucleonFcn.h 331 2013-05-13 13:39:22Z bartos $
 *
 * @file
 * @brief	Header for FCN function of nucleon FF.
 */

#ifndef NucleonFcn_H_
#define NucleonFcn_H_

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

#endif //NucleonFcn_H_
