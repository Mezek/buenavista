/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/FunctionGauss.h $
 * $Id$
 *
 * @file
 * @brief	Header for Gauss function.
 */

#ifndef _FunctionGauss_H_
#define _FunctionGauss_H_

namespace ROOT {

	namespace Minuit2 {

class FunctionGauss {

  private:

	double fMean;
	double fSigma;
	double fConstant;

  public:
  
	FunctionGauss(double mean, double sig, double constant) : 
		fMean(mean), fSigma(sig), fConstant(constant) {}

	double m() const {return fMean;}
	double s() const {return fSigma;}
	double c() const {return fConstant;}

	double operator()(double x) const {
    
		return c()*exp(-0.5*(x-m())*(x-m())/(s()*s()))/(sqrt(2.*TMath::Pi())*s());

	}
		
	~FunctionGauss() {}  
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // _FunctionGauss_H_
