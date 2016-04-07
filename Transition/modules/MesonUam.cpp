/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonUam.cpp $
 * $Id$
 *
 * @file
 * @brief	Model of U&A transition form factors: Eta, EtaPrime, PiZero.
 */

#include "MesonUam.h"

namespace ROOT {

	namespace Minuit2 {

FFactorT::FFactorT ()
{ 
	std::cout << ">> FFactor without defined size of parameters." << std::endl;

}

FFactorT::FFactorT ( int mp )
{
	modelPar = 5;
	a.resize(modelPar);
	v.resize(modelPar);
	allMesons = 3;
	switch (mp) {
		case 1: 
			FF.name = "Eta";
			FF.nor = 0.0348/massPi;
			FF.mesons = allMesons;
			break;
		case 2:
			FF.name = "EtaPrime";
			FF.nor = 0.0469/0.14;
			FF.mesons = allMesons;
			break;
		case 3:
			FF.name = "PiZero";
			FF.nor = 0.0352/0.14;
			FF.mesons = allMesons;
			break;
		default:
			FF.name = "Eta";
			FF.nor = 0.0348/0.14;
			FF.mesons = allMesons;
	}

	// Masses and widths
	mS[0] = massOm;
	mS[1] = massPhi;
	mS[2] = massOm1P;

	wS[0] = widthOm;
	wS[1] = widthPhi;
	wS[2] = widthOm1P;

	mV[0] = massRho;
	mV[1] = massRho1P;
	mV[2] = massRho2P;

	wV[0] = widthRho;
	wV[1] = widthRho1P;
	wV[2] = widthRho2P;
	
	for (int i = 0; i < 3; ++i) {
		mwS2[i] = (mS[i]-kI*wS[i]/2.)*(mS[i]-kI*wS[i]/2.);
		mwV2[i] = (mV[i]-kI*wV[i]/2.)*(mV[i]-kI*wV[i]/2.);
	}

	handSome = false;
	numberOfParameters = modelPar;
	cov.ResizeTo(modelPar, modelPar);
}

TComplex eL (const TComplex &a, const TComplex &b, const TComplex &c, const TComplex &sc, const double sign)
{
	TComplex f;
	if (sign == +1.) { f = (b-c)/(a-c)*(b-sc)/(a-sc)*(b-1./c)/(a-1./c)*(b-1./sc)/(a-1./sc); }
	if (sign == -1.) { f = (b-c)/(a-c)*(b-sc)/(a-sc)*(b+c)/(a+c)*(b+sc)/(a+sc); }
	return f;
}

/// Compare difference of two numbers with a small delta.

int howDiff ( double a, double b, double stdDiff )
{
	double p = fabs(a-b);
	//std::cout << ">> " << p << std::endl;
	if (p <= stdDiff) { return 0; } else { return 1; };
}

/// Load parameters from file

void FFactorT::LoadParameters (char* ds)
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

void FFactorT::SetParameters (const std::vector<double>& par)
{
	for (int i = 0; i < modelPar; ++i) {
		a[i].val = par[i];
		v[i].val = par[i];	
	}
}

/// Set i. parameter value

void FFactorT::SetParameter (const int i, const double value)
{
	a[i].val = value;
	v[i].val = value;	
}

/// Print parameters on std::cout

void FFactorT::PrintParameters ()
{
	std::cout << "\n>> Actual parameters of FFactor: " << FF.name << std::endl;
	for (int i = 0; i < modelPar; ++i) 
	{
		std::cout.width(3);
		std::cout << i+1 << ". ";
		std::cout.width(12);
		std::cout << this->v[i].name;
		std::cout.precision(5);
		std::cout.width(17);
		std::cout << this->v[i].val << "   +/- ";
		std::cout.width(17);
		std::cout << this->v[i].err << std::endl;
	}

	// Norms
	std::cout << "\n>> Norm:" << std::endl;
	std::cout << FF.name << ": " << FF.nor << std::endl;

	// Mesons under/over threshold
	std::cout << "\n>> Mesons placement to threshold:" << std::endl;
	std::cout.width(10);
	std::cout << "FF";
	std::cout.width(6);
	std::cout << "Under";
	std::cout.width(6);
	std::cout << "Over";
	std::cout.width(6);
	std::cout << "All";
	std::cout.width(6);
	std::cout << std::endl;

	for (int i = 0; i < 1; ++i) {
		int d = 0;
		int u = 0;
		double t;
		for (int j = 0; j < allMesons; ++j) {
			if ( i%2 == 0 ) {
				t = mS[j]*mS[j]-wS[j]*wS[j]/4.;
			} else {
				t = mV[j]*mV[j]-wV[j]*wV[j]/4.;
			}
			if ( t < a[i].val ) {
				d = d + 1;
			} else {
				u = u + 1;
			}
		}
		std::cout.width(10);
		std::cout << FF.name;
		std::cout.width(6);
		std::cout << d;
		std::cout.width(6);
		std::cout << u;
		std::cout.width(6);
		std::cout << u + d;
		std::cout << std::endl;
	}
}

/// Return value of i. parameter

double FFactorT::A ( int i )
{
	return a[i].val;
}

/// Return error of i. parameter

double FFactorT::E ( int i )
{
	return a[i].err;
}

/**
 * Basic W transformation:
 * sign = +1. : under threshold
 * sign = -1. : over threshold
 * It coressponds to change V(t) to -1./V(t)
 */

TComplex FFactorT::W (const TComplex &t, const TComplex &t0, const TComplex &tin, const double sign)
{
	TComplex cQ,cQin,cK,cZ,cRes;

	cQ = cQ.Sqrt((t-t0)/t0);
	cQin = cQin.Sqrt((tin-t0)/t0);

	//	bad!!!	if (cQ.Re() > tin.Re())
	if (t.Re() > tin.Re()) {
		cQ = cQ + 0.000001*kI;
	}

	cK = cK.Sqrt(cQin+cQ);
	cZ = cZ.Sqrt(cQin-cQ);
	cRes = kI*(cK-sign*cZ)/(cK+sign*cZ);
	return cRes;
}

/** Mesons:
    0 = Om
	1 = Phi
	2 = Om1P
*/

TComplex FFactorT::ScalarP (TComplex t)
{
	TComplex v = this->W(t,t0t,a[0].val,1.);
	TComplex vN = this->W(k0,t0t,a[0].val,1.);
	
	double sign, mt;
	for (int i = 0; i < FF.mesons; i++) {
		mt = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if ( mt < a[0].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = this->W(mwS2[i],t0t,a[0].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		//std::cout << sign << " ";
	}
	//std::cout << std::endl;

	TComplex norm, normA, suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,2);
	suma = normA*(0.5*FF.nor*mul[2] + (mul[0] - mul[2])*a[2].val +
		(mul[1] - mul[2])*a[3].val);
	return suma;
}

/** Mesons:
	0 = Rho
	1 = Rho1P
*/

TComplex FFactorT::VectorP (TComplex t)
{
	TComplex v = this->W(t,t0t,a[1].val,1.);
	TComplex vN = this->W(k0,t0t,a[1].val,1.);

	double sign, mt;
	for (int i = 0; i < FF.mesons; i++) {
		mt = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if ( mt < a[1].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = this->W(mwV2[i],t0t,a[1].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		//std::cout << sign << " ";
	}
	//std::cout << std::endl;	

	TComplex norm, normA, suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,2);	
	suma = normA*(0.5*FF.nor*mul[1] + (mul[0] - mul[1])*a[4].val);
	return suma;
}

TComplex FFactorT::FFVal ( TComplex t )
{
	TComplex val = this->ScalarP(t) + this->VectorP(t); 
	return val;
}

double FFactorT::FFAbsVal ( TComplex t )
{
	TComplex val = this->FFVal(t).Rho();
	return val;
}

	}  // namespace Minuit2

}  // namespace ROOT
