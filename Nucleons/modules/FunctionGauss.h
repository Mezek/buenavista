/**
 * $Date: 2013-05-13 15:39:22 +0200 (Mon, 13 May 2013) $
 * $Revision: 331 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/FunctionGauss.h $
 * $Id: FunctionGauss.h 331 2013-05-13 13:39:22Z bartos $
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
