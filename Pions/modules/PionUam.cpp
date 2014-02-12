/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Arbitrary UA form factor of pion + Omega-Rho interference.
 */

#include "PionUam.h"

namespace ROOT {

	namespace Minuit2 {

FFactor::FFactor ( std::size_t size ): a(size), v(size) 
{		
	t(1.,0.00001);
	//std::cout << "FF is empty! Value of t: " << t << std::endl;


	handSome = false;
	modelPar = a.size();
	numberOfParameters = modelPar;
	cov.ResizeTo(modelPar,modelPar);
}

/// Load parameters from file

void FFactor::LoadParameters (char* ds)
{
	ifstream myDataFile (ds);
	int k = 0;
	int N = 0;
	std::string line;
	if (myDataFile.is_open()) {
		while (myDataFile.peek() != EOF) {
			getline (myDataFile,line);
			std::stringstream stream (line);
			std::string field;
			int p = 0;
			std::string S = "";
			double A = 0.0;
			double B = 0.0;
			double C = 0.0;
			while (getline(stream,field,',')) {
				std::stringstream ss (field);
				if ( p == 0) { ss >> N; }
				if ( p == 1) { ss >> a[k].name; }
				if ( p == 2) { ss >> a[k].val; }
				if ( p == 3) { ss >> a[k].err; }
				if ( p == 4) { ss >> a[k].down; }
				if ( p == 5) { ss >> a[k].up; }
				++p;
			}
			//std::cout << "N: " << N << std::endl;
			
			this->v[k].name = a[k].name;
			this->v[k].val = a[k].val;
			this->v[k].err = a[k].err;
			this->v[k].down = a[k].down;
			this->v[k].up = a[k].up;
			//std::cout << N << ".: " << v[k].name << " : " << v[k].val << " : " << v[k].err << " : " << v[k].down << " : " << v[k].down << std::endl;
			++k;
		}
		myDataFile.close();
		if (N != modelPar) {
			std::cout << "> LoadParameters: Error!" << std::endl;
			std::cout << "> Number of parameters: " << N << " in '" << ds << "' differs from model: " << modelPar << "!" << std::endl;
		}
	}
	else std::cerr << ">> LoadParameters: Error! Unable to open parametric file: '" << ds << "'!" << std::endl;

}

/// Set parameters from other source than file

void FFactor::SetParameters (const std::vector<double>& par)
{
	for (int i = 0; i < modelPar; ++i) {
		a[i].val = par[i];
	}
}

/// Set i. parameter value

void FFactor::SetParameter (const int i, const double value)
{
	a[i].val = value;
}

/// Print parameters on std::cout

void FFactor::PrintParameters ()
{
	std::cout << "\n>> Actual parameters of FFactor:" << std::endl;
	for (int i = 0; i < modelPar; ++i) 
	{
		std::cout.width(3);
		std::cout << i+1 << ". ";
		std::cout.width(12);
		std::cout << this->a[i].name;
		std::cout.precision(10);
		std::cout.width(17);
		std::cout << this->a[i].val << "   +/- ";
		std::cout.width(17);
		std::cout << this->a[i].err << std::endl;
	}
}

/// Return value of i. parameter

double FFactor::A ( int i )
{
	return a[i].val;
}

/// Return error of i. parameter

double FFactor::E ( int i )
{
	return a[i].err;
}

/// Load covariance matrix from file
	
void FFactor::LoadCovMatrix ( char* d )
{
	ifstream myDataFile (d);
	if (myDataFile.is_open()) {
		int num = 0;
		char firstChar;
		std::string line;
		double dA;
		//cov.resize(modelPar);
		while (myDataFile.peek() != EOF) {
			firstChar = myDataFile.peek();
			if ( (firstChar == '%') || (firstChar == '#') || (firstChar == '*') || (firstChar == '/') ) {
				getline (myDataFile,line);
			}
			else {
				for (int i = 0; i < modelPar; ++i) {
					myDataFile >> dA;
					cov(num,i) = dA;
				}
				getline (myDataFile,line);
				++num;
			}
		}
		myDataFile.close();
	}
	else std::cerr << "> LoadCovMatrix: Error! Unable to open covariance matrix file '" << d << "'" << std::endl;
}

/**
 * Basic W transformation:
 * sign = +1. : under threshold
 * sign = -1. : over threshold
 * It coressponds to change V(t) to -1./V(t)
 */

TComplex FFactor::W (const TComplex &t, const double sign) {
	TComplex cA,ct0v,cQ,cQin,cK,cZ,cRes;
	cA(a[0],0);
	ct0v(t0v,0);

	cQ = cQ.Sqrt((t-ct0v)/ct0v);

	if (t.Re() > a[0]) { // bad: cQ.Re() > a[0]!
		cQ = cQ + 0.000001*kI;
		//std::cout << "cQ: " << a[0] << " " << t << " " << cQ << std::endl;
	}

	cQin = cQin.Sqrt((cA-ct0v)/ct0v);
	cK = cK.Sqrt(cQin+cQ);
	cZ = cZ.Sqrt(cQin-cQ);
	cRes = kI*(cK-sign*cZ)/(cK+sign*cZ);
	return cRes;
}

TComplex FFactor::Value (TComplex t) {
	TComplex cWn,cW;
	TComplex cQrho,cW1,cW1c;
	TComplex cQrho1,cW2,cW2c;
	TComplex cQrho2,cW3,cW3c;
	TComplex cNrho,cNrho1,cNrho2;
	TComplex cDrho,cDrho1,cDrho2;
	TComplex cfw,cf1w,cf2w,cf3w,cf4w;
	TComplex pkv0,pkv1,pkv2,pkv3,pkvc;
	TComplex c2Sq,c3Sq,cNW1,cNW2,pkv2a,pkv3a;
	TComplex cNorm,cMult,cFpi,cMomr,cEp,cSuma;
	double fazi;

	cWn = FFactor::W(k0,1.);
	cW = FFactor::W(t,1.);

	cQrho = (a[1]-kI*a[2]/2.)*(a[1]-kI*a[2]/2.);
	cW1 = FFactor::W(cQrho,1.);
	cW1c = cW1c.Conjugate(cW1);

	cQrho1 = (a[3]-kI*a[4]/2.)*(a[3]-kI*a[4]/2.);
	cW2 = FFactor::W(cQrho1,-1.);
	cW2c = cW2c.Conjugate(cW2);
	
	cQrho2 = (a[5]-kI*a[6]/2.)*(a[5]-kI*a[6]/2.);
	cW3 = FFactor::W(cQrho2,-1.);
	cW3c = cW3c.Conjugate(cW3);

	cNrho = (cWn-cW1)*(cWn-cW1c)*(cWn-1./cW1)*(cWn-1./cW1c);
	cDrho = (cW-cW1)*(cW-cW1c)*(cW-1./cW1)*(cW-1./cW1c);
	cNrho1 = (cWn-cW2)*(cWn-cW2c)*(cWn+cW2)*(cWn+cW2c);
	cDrho1 = (cW-cW2)*(cW-cW2c)*(cW+cW2)*(cW+cW2c);
	cNrho2 = (cWn-cW3)*(cWn-cW3c)*(cWn+cW3)*(cWn+cW3c);
	cDrho2 = (cW-cW3)*(cW-cW3c)*(cW+cW3)*(cW+cW3c);
	
	cfw = cNrho/cDrho;
	cf1w = cNrho1/cDrho1;
	cf2w = cNrho2/cDrho2;
	cf3w = cNrho/cNrho2;
	cf4w = cNrho1/cNrho2;

	c2Sq = (cW2*cW2c);
	c3Sq = (cW3*cW3c);
	cNW1 = cNrho1/c2Sq/c2Sq;
	cNW2 = cNrho2/c3Sq/c3Sq;

	pkv0 = pkv0.Power((cW3*cW3c)/(cW2*cW2c),2);
	pkv1 = 1.+(a[9]*a[10])/(a[9]-a[10])*2.*cW1.Re()*(1.+1./(cW1*cW1c));
	pkv2 = (1.-(1.-pkv1*c3Sq*c3Sq*cf3w)*a[7])/(1.-pkv0*cf4w);
	pkv3 = (-1.)*pkv1*c3Sq*c3Sq*cf3w*a[7] - pkv0*cf4w*pkv2;
	pkvc = (a[7] + pkv2 + pkv3);

	// Alternative constants
	pkv3a = (cNW1-(cNW1-pkv1*cNrho)*a[7])/(cNW1-cNW2);
	pkv2a = (-cNW2+(cNW2-pkv1*cNrho)*a[7])/(cNW1-cNW2);
	
	cNorm = cNorm.Power((1.-cW*cW)/(1.-cWn*cWn),2);
	cMult = (cW-a[9])*(cWn-a[10])/(cWn-a[9])/(cW-a[10]);
	cFpi = cNorm*(cfw*a[7]+cf1w*pkv2+cf2w*pkv3)*cMult;
	fazi = atan2((a[1]*a[2]),(a[1]*a[1]-massOm*massOm));
	cMomr = massOm*massOm - t - kI*massOm*widthOm;
	cEp = cEp.Exp(kI*fazi);
	
	cSuma = cFpi + a[8]*cEp*massOm*massOm/cMomr;

	return cSuma;
}

double FFactor::ValueSquared (TComplex t) {
	TComplex vSq = FFactor::value(t);
	return vSq.Re()*vSq.Re()+vSq.Im()*vSq.Im();
}

	}  // namespace Minuit2

}  // namespace ROOT
#endif // _PionUam_H_
