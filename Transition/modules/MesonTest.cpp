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

#include "MesonTest.h"

namespace ROOT {

	namespace Minuit2 {
		
FFactor::FFactor ( std::size_t size ): a(size), v(size)
{ 
	t(1.,0.00001);
	std::cout << ">> FFactor is empty! Value of t: " << t << std::endl;

	// izo :: name, nor, mesons

	// Names Dirac, Pauli / Sachs
	FF[0].name = "F1s";
	FF[1].name = "F1v";
	FF[2].name = "F2s";
	FF[3].name = "F2v";

	FFtype[0] = "GEP";
	FFtype[1] = "GMP";
	FFtype[2] = "GEN";
	FFtype[3] = "GMN";
		
	// Norms
	FF[0].nor = 0.5;
	FF[1].nor = 0.5;
	FF[2].nor = 0.5*(ammP+ammN);
	FF[3].nor = 0.5*(ammP-ammN);

	// All mesons in model
	FF[0].mesons = 6;
	FF[1].mesons = 3;
	FF[2].mesons = 6;
	FF[3].mesons = 3;
	allMesons = 18;

	// Masses and widths
	mS[0] = massOm;
	mS[1] = massPhi;
	mS[2] = massOm1P;
	mS[3] = massPhi1P;
	mS[4] = massOm2P;
	mS[5] = massPhi2P;

	wS[0] = widthOm;
	wS[1] = widthPhi;
	wS[2] = widthOm1P;
	wS[3] = widthPhi1P;
	wS[4] = widthOm2P;
	wS[5] = widthPhi2P;
	
	mV[0] = massRho;
	mV[1] = massRho1P;
	mV[2] = massRho2P;
	mV[3] = massRho3P;
	mV[4] = 0.;
	mV[5] = 0.;

	wV[0] = widthRho;
	wV[1] = widthRho1P;
	wV[2] = widthRho2P;
	wV[3] = widthRho3P;
	wV[4] = 0.;
	wV[5] = 0.;
	
	for (int i = 0; i < 6; ++i) {
		mwS2[i] = (mS[i]-kI*wS[i]/2.)*(mS[i]-kI*wS[i]/2.);
		mwV2[i] = (mV[i]-kI*wV[i]/2.)*(mV[i]-kI*wV[i]/2.);
	}

	handSome = false;
	modelPar = a.size();
	numberOfParameters = modelPar;
	cov.ResizeTo(modelPar, modelPar);
	expressPar = allMesons - modelPar + 4;
}

TComplex eL (const TComplex &a, const TComplex &b, const TComplex &c, const TComplex &sc, const double sign)
{
	TComplex f;
	if (sign == +1.) { f = (b-c)/(a-c)*(b-sc)/(a-sc)*(b-1./c)/(a-1./c)*(b-1./sc)/(a-1./sc); }
	if (sign == -1.) { f = (b-c)/(a-c)*(b-sc)/(a-sc)*(b+c)/(a+c)*(b+sc)/(a+sc); }
	return f;
}

TComplex sI (const TComplex &b, const TComplex &c, const TComplex &sc, const double sign)
{
	TComplex f;
	if (sign == +1.) { f = -1.*(b-c)*(b-sc)/(c-1./c)*(b-1./c)*(b-1./sc)/(sc-1./sc); }
	if (sign == -1.) { f = -1.*(b-c)*(b-sc)/(c-1./c)*(b+c)*(b+sc)/(sc-1./sc); }
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
		v[i].val = par[i];	
	}
}

/// Set i. parameter value

void FFactor::SetParameter (const int i, const double value)
{
	a[i].val = value;
	v[i].val = value;	
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
		std::cout << this->v[i].name;
		std::cout.precision(5);
		std::cout.width(17);
		std::cout << this->v[i].val << "   +/- ";
		std::cout.width(17);
		std::cout << this->v[i].err << std::endl;
	}

	// Norms
	std::cout << "\n>> Norms for the components:" << std::endl;
	for (int i = 0; i < 4; ++i) {
		std::cout << FF[i].name << ": nor[" << i << "] = " << FF[i].nor << std::endl;
	}

	// Mesons under/over threshold
	std::cout << "\n>> Mesons placement to threshold:" << std::endl;
	std::cout << "       Under  Over  All" << std::endl;

	for (int i = 0; i < 4; ++i) {
		int d = 0;
		int u = 0;
		double t;
		for (int j = 0; j < FF[i].mesons; ++j) {
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
		std::cout.width(4);
		std::cout << FF[i].name;
		std::cout.width(6);
		std::cout << d;
		std::cout.width(6);
		std::cout << u;
		std::cout.width(6);
		std::cout << u + d;
		std::cout << std::endl;
	}
}



	}  // namespace Minuit2

}  // namespace ROOT
