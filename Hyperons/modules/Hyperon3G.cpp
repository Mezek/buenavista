/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/modules/Nucleon3GH.cpp $
 * $Id$
 *
 * @file
 * @brief	U&A form factor for nucleons &mdash; 3 vector meson resonances, v. hyperon.
 */

#include "NucleonUamH.h"

/**
 * Namespace for some useful functions.
 */

namespace shortCut {



/// Signum of given number.

std::string signum ( double n )
{
	std::string s;
	if ( (n > 0.) || (n=0.) ) {
		s = "+";
	} else {
		s = "-";
	}
	return s;
}

}  // namespace shortCut


namespace ROOT {

	namespace Minuit2 {

	// Explicitly define number of parameters for used FF
	const int globalFF = 13;

FFactor::FFactor (): a(globalFF), b(globalFF), v(globalFF), particleType(1)
{
	FFactor::SetParticle(particleType);
}

FFactor::FFactor ( int myParticle ): a(globalFF), b(globalFF), v(globalFF), particleType(myParticle)
{
	FFactor::SetParticle(particleType);
}

void FFactor::SetParticle ( int myParticle )
{
	//std::cout << "\n>> FFactor has type: " << myParticle << std::endl;
	particleType = myParticle;
	
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

	// izo :: name, nor, mesons

	// Names
	FF[0].name = "F1s";
	FF[1].name = "F1v";
	FF[2].name = "F2s";
	FF[3].name = "F2v";

	// Norms
	while ( HyperonsCheckRange(particleType) != 1 ) {
		particleType = HyperonsGetType();
		std::cout << ">> FFactor set for particle type: " << particleType << std::endl;
	}

	switch (particleType) {
		case 1:
			FF[0].nor = 0.5;
			FF[1].nor = 0.5;
			FF[2].nor = 0.5*(muP+muN-1.);
			FF[3].nor = 0.5*(muP-muN-1.);
			massH1 = massP;
			massH2 = massN;
			particleName1 = "p";
			particleName2 = "n";
			break;
		case 2:
			FF[0].nor = 0.;
			FF[1].nor = 0.;
			FF[2].nor = muL;
			FF[3].nor = 0.;
			massH1 = massL;
			massH2 = 0.;
			particleName1 = "Lambda";
			particleName2 = "-";
			break;
		case 3:
			FF[0].nor = 0.;
			FF[1].nor = 1.;
			FF[2].nor = 0.5*(muSp+muSm);
			FF[3].nor = 0.5*(muSp-muSm-2.);
			massH1 = massSp;
			massH2 = massSm;
			particleName1 = "Sigma+";
			particleName2 = "Sigma-";
			break;
		case 4:
			FF[0].nor = 0.;
			FF[1].nor = 0.;
			FF[2].nor = 0.5*(muSp+muSm);
			FF[3].nor = 0.;
			massH1 = massS0;
			massH2 = 0.;
			particleName1 = "Sigma0";
			particleName2 = "-";
			break;			
		case 5:
			FF[0].nor = -0.5;
			FF[1].nor = 0.5;
			FF[2].nor = 0.5*(muX0+muXm+1.);
			FF[3].nor = 0.5*(muX0-muXm-1.);
			massH1 = massX0;
			massH2 = massXm;
			particleName1 = "Xi0";
			particleName2 = "Xim";
			break;
		default:
			std::cout << ">> Error: No particle type defined!" << std::endl;
	}

	// All mesons in model
	FF[0].mesons = 6;
	FF[1].mesons = 3;
	FF[2].mesons = 6;
	FF[3].mesons = 3;
	allMesons = 18;

	handSome = false;
	modelPar = globalFF;
	numberOfParameters = modelPar;
	cov.ResizeTo(modelPar, modelPar);
	expressPar = 10;

}

/// Get particle type

int FFactor::GetParticle( void )
{
	return particleType;
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

	this->FixParameters();
	// debug: if ( particleType == 1 ) { particleType = 5; }	
	if ( particleType != 1 ) { this->TransformSU3(); }
}

void FFactor::FixParameters ()
{
	for (int i = 0; i < modelPar; ++i) {
		v[i].name = a[i].name;
		v[i].val = a[i].val;
		v[i].err = a[i].err;
		v[i].down = a[i].down;
		v[i].up = a[i].up;
	}
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
	this->ExpressedParameters();
	std::cout << std::endl;
	for (int i = 0; i < expressPar; i++) {
		std::cout.width(3);
		std::cout << i+1 << ".";
		std::cout.width(10);
		std::cout << b[i].name;
		std::cout.width(13);
		std::cout << b[i].val << std::endl;
	}
}

void FFactor::ExpressedParameters ()
{
	std::cout << "\n>> Calculated expressed parameters:" << std::endl;
	std::vector<hod> z(allMesons);
	double sign, mt;
	TComplex vN;

	// F1s
	vN = FFactor::W(k0,t0s,a[0].val,1.);
	for (int i = 0; i < FF[0].mesons; i++) {
		mt = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if ( mt < a[0].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = FFactor::W(mwS2[i],t0s,a[0].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}
	b[0].name = "f_Om2";
	b[0].val = sub[5]/(sub[5]-sub[4])*0.5 
		- (sub[5]-sub[0])/(sub[5]-sub[4])*a[4].val
		- (sub[5]-sub[1])/(sub[5]-sub[4])*a[5].val
		- (sub[5]-sub[2])/(sub[5]-sub[4])*a[6].val
		- (sub[5]-sub[3])/(sub[5]-sub[4])*a[7].val;

	b[1].name = "f_Ph2";
	b[1].val = - sub[4]/(sub[5]-sub[4])*0.5 
		+ (sub[4]-sub[0])/(sub[5]-sub[4])*a[4].val
		+ (sub[4]-sub[1])/(sub[5]-sub[4])*a[5].val
		+ (sub[4]-sub[2])/(sub[5]-sub[4])*a[6].val
		+ (sub[4]-sub[3])/(sub[5]-sub[4])*a[7].val;

	// F2s
	vN = FFactor::W(k0,t0s,a[2].val,1.);
	for (int i = 0; i < FF[2].mesons; i++) {
		mt = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if ( mt < a[2].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = FFactor::W(mwS2[i],t0s,a[2].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}

	b[2].name = "f_Om1_T";
	b[2].val = sub[5]/(sub[5]-sub[2])*sub[4]/(sub[4]-sub[2])*0.5*(muP+muN-1.)
		- (sub[5]-sub[0])/(sub[5]-sub[2])*(sub[4]-sub[0])/(sub[4]-sub[2])*a[9].val
		- (sub[5]-sub[1])/(sub[5]-sub[2])*(sub[4]-sub[1])/(sub[4]-sub[2])*a[10].val
		- (sub[5]-sub[3])/(sub[5]-sub[2])*(sub[4]-sub[3])/(sub[4]-sub[2])*a[11].val;
	b[3].name = "f_Om2_T";
	b[3].val = - sub[5]/(sub[5]-sub[4])*sub[2]/(sub[4]-sub[2])*0.5*(muP+muN-1.)
		+ (sub[5]-sub[0])/(sub[5]-sub[4])*(sub[2]-sub[0])/(sub[4]-sub[2])*a[9].val
		+ (sub[5]-sub[1])/(sub[5]-sub[4])*(sub[2]-sub[1])/(sub[4]-sub[2])*a[10].val
		+ (sub[5]-sub[3])/(sub[5]-sub[4])*(sub[2]-sub[3])/(sub[4]-sub[2])*a[11].val;
	b[4].name = "f_Ph2_T";
	b[4].val = sub[2]/(sub[5]-sub[2])*sub[4]/(sub[5]-sub[4])*0.5*(muP+muN-1.)
		- (sub[2]-sub[0])/(sub[5]-sub[2])*(sub[4]-sub[0])/(sub[5]-sub[4])*a[9].val
		- (sub[2]-sub[1])/(sub[5]-sub[2])*(sub[4]-sub[1])/(sub[5]-sub[4])*a[10].val
		- (sub[2]-sub[3])/(sub[5]-sub[2])*(sub[4]-sub[3])/(sub[5]-sub[4])*a[11].val;

/*	// choice Bartos
	b[2].name = "f_Om1_T";
	b[2].val = sub[5]/(sub[5]-sub[3])*sub[4]/(sub[4]-sub[3])*0.5*(muP+muN-1.) 
		- (sub[5]-sub[0])/(sub[5]-sub[3])*(sub[4]-sub[0])/(sub[4]-sub[3])*a[9].val
		- (sub[5]-sub[1])/(sub[5]-sub[3])*(sub[4]-sub[1])/(sub[4]-sub[3])*a[10].val
		- (sub[5]-sub[2])/(sub[5]-sub[3])*(sub[4]-sub[2])/(sub[4]-sub[3])*a[11].val;
	b[3].name = "f_Om2_T";
	b[3].val = - sub[5]/(sub[5]-sub[4])*sub[3]/(sub[4]-sub[3])*0.5*(muP+muN-1.)
		+ (sub[5]-sub[0])/(sub[5]-sub[4])*(sub[3]-sub[0])/(sub[4]-sub[3])*a[9].val
		+ (sub[5]-sub[1])/(sub[5]-sub[4])*(sub[3]-sub[1])/(sub[4]-sub[3])*a[10].val
		+ (sub[5]-sub[2])/(sub[5]-sub[4])*(sub[3]-sub[2])/(sub[4]-sub[3])*a[11].val;
	b[4].name = "f_Ph2_T";
	b[4].val = sub[3]/(sub[5]-sub[3])*sub[4]/(sub[5]-sub[4])*0.5*(muP+muN-1.)
		- (sub[3]-sub[0])/(sub[5]-sub[3])*(sub[4]-sub[0])/(sub[5]-sub[4])*a[9].val
		- (sub[3]-sub[1])/(sub[5]-sub[3])*(sub[4]-sub[1])/(sub[5]-sub[4])*a[10].val
		- (sub[3]-sub[2])/(sub[5]-sub[3])*(sub[4]-sub[2])/(sub[5]-sub[4])*a[11].val;
*/
	
	// F1v
	vN = this->W(k0,t0v,a[1].val,1.);
	for (int i = 0; i < FF[1].mesons; i++) {
		mt = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if ( mt < a[1].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = FFactor::W(mwV2[i],t0v,a[1].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}
	b[5].name = "f_Rh1";
	b[5].val = sub[2]/(sub[2]-sub[1])*0.5 
		- (sub[2]-sub[0])/(sub[2]-sub[1])*a[8].val;

	b[6].name = "f_Rh2";
	b[6].val = sub[0]/(sub[0]-sub[2])*0.5 
		- (sub[0]-sub[1])/(sub[0]-sub[2])*b[5].val;

	// F2v
	vN = FFactor::W(k0,t0v,a[3].val,1.);
	for (int i = 0; i < FF[3].mesons; i++) {
		mt = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if ( mt < a[3].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = FFactor::W(mwV2[i],t0v,a[3].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}
	b[7].name = "f_Rh_T";
	b[7].val = sub[1]/(sub[1]-sub[0])*sub[2]/(sub[2]-sub[0])*0.5*(muP-muN-1.);
	b[8].name = "f_Rh1_T";
	b[8].val = sub[0]/(sub[0]-sub[1])*sub[2]/(sub[2]-sub[1])*0.5*(muP-muN-1.);
	b[9].name = "f_Rh2_T";
	b[9].val = sub[0]/(sub[0]-sub[2])*sub[1]/(sub[1]-sub[2])*0.5*(muP-muN-1.);
}

double Duo ( double &a, double &b )
{
	return a/(a-b);
}

double Trio ( double &a, double &b, double &c )
{
	return (a-b)/(a-c);
}

void FFactor::TransformSU3 ()
{
	// Calculate SU(3) relations
	double th[3];
	th[0] = 43.8*TMath::DegToRad();
	th[1] = 50.3*TMath::DegToRad();
	th[2] = 50.3*TMath::DegToRad();
	double k1, k2, k3, k4;

	double fuO[3], fuP[3], fuR[3]; // universal coupling constants
	fuO[0] = 17.0620;
	fuP[0] = -13.4428;
	fuR[0] = 4.9582;
	fuO[1] = 47.6022;
	fuP[1] = -33.6598; //37.0510;
	fuR[1] = 13.6491;
	fuO[2] = 48.3651;
	fuP[2] = -34.1993;
	fuR[2] = 22.5275;

	this->ExpressedParameters();

	// first index for i=(1),(2), second index for family j=1,2,3
	double fO[2][3], fP[2][3], fR[2][3];
	fO[0][0] =  a[4].val*fuO[0];
	fP[0][0] =  a[5].val*fuP[0];
	fR[0][0] =  a[8].val*fuR[0];
	fO[1][0] =  a[9].val*fuO[0];
	fP[1][0] = a[10].val*fuP[0];
	fR[1][0] =  b[7].val*fuR[0];

	fO[0][1] =  a[6].val*fuO[1];
	fP[0][1] =  a[7].val*fuP[1];
	fR[0][1] =  b[5].val*fuR[1];
	fO[1][1] =  b[2].val*fuO[1];
	fP[1][1] = a[11].val*fuP[1];
	fR[1][1] =  b[8].val*fuR[1];

	fO[0][2] = fuO[2];
	fP[0][2] = fuP[2];
	fR[0][2] = fuR[2];
	fO[1][2] = fuO[2];
	fP[1][2] = fuP[2];
	fR[1][2] = fuR[2];

	// Inverse equations
	double fD[2][3], fF[2][3], fS[2][3];
	double fOp[2][3], fPp[2][3], fRp[2][3];
	std::cout << "Generation " << "  fD  " << "     fF  " << "       fS" << std::endl;
	for (int j = 0; j < 3; ++j) {
		for (int i = 0; i < 2; ++i) {
			fD[i][j] = TMath::Sqrt(3.)/2.*(fO[i][j]*sin(th[j]) - fP[i][j]*cos(th[j])) + 1.5*fR[i][j];
			fF[i][j] = TMath::Sqrt(3.)/2.*(-fO[i][j]*sin(th[j]) + fP[i][j]*cos(th[j])) + 0.5*fR[i][j];
			fS[i][j] = TMath::Sqrt(2.)*(fO[i][j]*cos(th[j]) + fP[i][j]*sin(th[j]));
			std::cout.precision(6);
			std::cout.width(7);
			std::cout << j << "   " << fD[i][j] << "   " << fF[i][j] << "   " << fS[i][j] << std::endl;

			switch (particleType) {
				case 1: // just for completness
					fOp[i][j] = fO[i][j];
					fPp[i][j] = fP[i][j];
					fRp[i][j] = fR[i][j];
					break;
				case 2:
					fOp[i][j] = k1*fS[i][j] + k2*fD[i][j];
					fPp[i][j] = k3*fS[i][j] - k4*fD[i][j];
					fRp[i][j] = 0.;
					break;
				case 3:
					fOp[i][j] = k1*fS[i][j] - k3*fD[i][j];
					fPp[i][j] = k3*fS[i][j] + k4*fD[i][j];
					fRp[i][j] = -fF[i][j];
					break;
				case 4:
					fOp[i][j] = k1*fS[i][j] - k3*fD[i][j];
					fPp[i][j] = k3*fS[i][j] + k4*fD[i][j];
					fRp[i][j] = -fF[i][j];
					break;			
				case 5:
					fOp[i][j] = k1*fS[i][j] + 3./2.*k2*fF[i][j] + k2/2.*fD[i][j];
					fPp[i][j] = k3*fS[i][j] - 3./2.*k4*fF[i][j] - k4/2.*fD[i][j];
					fRp[i][j] = (fD[i][j] - fF[i][j])/2.;
					break;
				default:
					std::cout << ">> Error: No particle type defined!" << std::endl;
			}
		}
	}

	a[4].val = fOp[0][0]/fuO[0];
	a[5].val = fPp[0][0]/fuP[0];
	a[6].val = fOp[0][1]/fuO[1];
	a[7].val = fPp[0][1]/fuP[1];
	a[8].val = fRp[0][0]/fuR[0];
	a[9].val = fOp[1][0]/fuO[0];
	a[10].val = fPp[1][0]/fuP[0];
	a[11].val = fPp[1][1]/fuP[1];
	this->FixParameters();	
	std::cout << ">> SU(3) transformation of coupling constants... done." << std::endl;

	std::cout.precision(6);
	std::cout << "1 f_Om : " << fO[0][0] << std::endl;
	std::cout << "1 f_Ph : " << fP[0][0] << std::endl;
	std::cout << "1 f_Rh : " << fR[0][0] << std::endl;
	std::cout << "T f_Om : " << fO[1][0] << std::endl;
	std::cout << "T f_Ph : " << fP[1][0] << std::endl;
	std::cout << "T f_Rh : " << fR[1][0] << std::endl;

	std::cout << "1 f_Om1: " << fO[0][1] << std::endl;
	std::cout << "1 f_Ph1: " << fP[0][1] << std::endl;
	std::cout << "1 f_Rh1: " << fR[0][1] << std::endl;
	std::cout << "T f_Om1: " << fO[1][1] << std::endl;
	std::cout << "T f_Ph1: " << fP[1][1] << std::endl;
	std::cout << "T f_Rh1: " << fR[1][1] << std::endl;
	this->PrintParameters();
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

TComplex FFactor::W (const TComplex &t, const TComplex &t0, const TComplex &tin, const double sign)
{
	TComplex cQ,cQin,cK,cZ,cRes;

	cQ = cQ.Sqrt((t-t0)/t0);
	cQin = cQin.Sqrt((tin-t0)/t0);

	//	bad!!!	if (cQ.Re() > tin.Re())
	if (t.Re() > tin.Re()) {
		cQ = cQ + 0.000001*kI;
		//std::cout << ">> t: "<< t.Re() << " " << tin.Re() << " " << cQ << std::endl;
	}

	cK = cK.Sqrt(cQin+cQ);
	cZ = cZ.Sqrt(cQin-cQ);
	cRes = kI*(cK-sign*cZ)/(cK+sign*cZ);
	return cRes;
}

/** 0 = Om, 1 = Phi, 2 = Om1P, 3 = Phi1P, 4 = Om2P, 5 = Phi2P */

TComplex FFactor::ScalarOne (TComplex t)
{
	TComplex v = FFactor::W(t,t0s,a[0].val,1.);
	TComplex vN = FFactor::W(k0,t0s,a[0].val,1.);
	
	double sign, mt;
	for (int i = 0; i < FF[0].mesons; i++) {
		mt = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if ( mt < a[0].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = FFactor::W(mwS2[i],t0s,a[0].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}
	//std::cout << ">> 0-2: "<< mul[0] << " " << mul[1] << " " << mul[2] << std::endl;
	//std::cout << ">> 3-6: "<< mul[3] << " " << mul[4] << " " << mul[5] <<  " " << mul[6] << std::endl;

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,4);
	suma = normA*(FF[0].nor*mul[4]*mul[5] +
	  (mul[0]*mul[5]*(sub[5]-sub[0])/(sub[5]-sub[4]) +
	   mul[0]*mul[4]*(sub[4]-sub[0])/(sub[4]-sub[5]) -
	   mul[4]*mul[5])*a[4].val +
	  (mul[1]*mul[5]*(sub[5]-sub[1])/(sub[5]-sub[4]) +
	   mul[1]*mul[4]*(sub[4]-sub[1])/(sub[4]-sub[5]) -
	   mul[4]*mul[5])*a[5].val +
	  (mul[2]*mul[5]*(sub[5]-sub[2])/(sub[5]-sub[4]) +
	   mul[2]*mul[4]*(sub[4]-sub[2])/(sub[4]-sub[5]) -
	   mul[4]*mul[5])*a[6].val +
	  (mul[3]*mul[5]*(sub[5]-sub[3])/(sub[5]-sub[4]) +
	   mul[3]*mul[4]*(sub[4]-sub[3])/(sub[4]-sub[5]) -
	   mul[4]*mul[5])*a[7].val);

	//std::cout << ">> ScalarOne: "<< norm << " " << FF[0].nor << std::endl;
	return suma;
}

/** 0 = Om, 1 = Phi, 2 = Om1P, 3 = Phi1P, 4 = Om2P, 5 = Phi2P */

TComplex FFactor::ScalarTwo (TComplex t)
{
	TComplex v = FFactor::W(t,t0s,a[2].val,1.);
	TComplex vN = FFactor::W(k0,t0s,a[2].val,1.);
	
	double sign, mt;
	for (int i = 0; i < FF[2].mesons; i++) {
		mt = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if ( mt < a[2].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = FFactor::W(mwS2[i],t0s,a[2].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,6);
	// choice Dubnicka
	suma = normA*(FF[2].nor*mul[2]*mul[4]*mul[5] +
	  (mul[0]*mul[4]*mul[5]*(sub[4]-sub[0])/(sub[4]-sub[2])*(sub[5]-sub[0])/(sub[5]-sub[2]) +
	   mul[0]*mul[2]*mul[5]*(sub[2]-sub[0])/(sub[2]-sub[4])*(sub[5]-sub[0])/(sub[5]-sub[4]) +
	   mul[0]*mul[2]*mul[4]*(sub[2]-sub[0])/(sub[2]-sub[5])*(sub[4]-sub[0])/(sub[4]-sub[5]) -
	   mul[2]*mul[4]*mul[5])*a[9].val +
	  (mul[1]*mul[4]*mul[5]*(sub[4]-sub[1])/(sub[4]-sub[2])*(sub[5]-sub[1])/(sub[5]-sub[2]) +
	   mul[1]*mul[2]*mul[5]*(sub[2]-sub[1])/(sub[2]-sub[4])*(sub[5]-sub[1])/(sub[5]-sub[4]) +
	   mul[1]*mul[2]*mul[4]*(sub[2]-sub[1])/(sub[2]-sub[5])*(sub[4]-sub[1])/(sub[4]-sub[5]) -
	   mul[2]*mul[4]*mul[5])*a[10].val +
	  (mul[3]*mul[4]*mul[5]*(sub[4]-sub[3])/(sub[4]-sub[2])*(sub[5]-sub[3])/(sub[5]-sub[2]) +
	   mul[3]*mul[2]*mul[5]*(sub[2]-sub[3])/(sub[2]-sub[4])*(sub[5]-sub[3])/(sub[5]-sub[4]) +
	   mul[3]*mul[2]*mul[4]*(sub[2]-sub[3])/(sub[2]-sub[5])*(sub[4]-sub[3])/(sub[4]-sub[5]) -
	   mul[2]*mul[4]*mul[5])*a[11].val);

/*	// choice Bartos
	suma = normA*(FF[2].nor*mul[3]*mul[4]*mul[5] +
	  (mul[0]*mul[4]*mul[5]*(sub[4]-sub[0])/(sub[4]-sub[3])*(sub[5]-sub[0])/(sub[5]-sub[3]) +
	   mul[0]*mul[3]*mul[5]*(sub[3]-sub[0])/(sub[3]-sub[4])*(sub[5]-sub[0])/(sub[5]-sub[4]) +
	   mul[0]*mul[3]*mul[4]*(sub[3]-sub[0])/(sub[3]-sub[5])*(sub[4]-sub[0])/(sub[4]-sub[5]) -
	   mul[3]*mul[4]*mul[5])*a[9].val +
	  (mul[1]*mul[4]*mul[5]*(sub[4]-sub[1])/(sub[4]-sub[3])*(sub[5]-sub[1])/(sub[5]-sub[3]) +
	   mul[1]*mul[3]*mul[5]*(sub[3]-sub[1])/(sub[3]-sub[4])*(sub[5]-sub[1])/(sub[5]-sub[4]) +
	   mul[1]*mul[3]*mul[4]*(sub[3]-sub[1])/(sub[3]-sub[5])*(sub[4]-sub[1])/(sub[4]-sub[5]) -
	   mul[3]*mul[4]*mul[5])*a[10].val +
	  (mul[2]*mul[4]*mul[5]*(sub[4]-sub[2])/(sub[4]-sub[3])*(sub[5]-sub[2])/(sub[5]-sub[3]) +
	   mul[2]*mul[3]*mul[5]*(sub[3]-sub[2])/(sub[3]-sub[4])*(sub[5]-sub[2])/(sub[5]-sub[4]) +
	   mul[2]*mul[3]*mul[4]*(sub[3]-sub[2])/(sub[3]-sub[5])*(sub[4]-sub[2])/(sub[4]-sub[5]) -
	   mul[3]*mul[4]*mul[5])*a[11].val);
*/
	return suma;
}

/** 0 = Rho, 1 = Rho1P, 2 = Rho2P, 3 = Rho3P */

TComplex FFactor::VectorOne (TComplex t)
{
	TComplex v = FFactor::W(t,t0v,a[1].val,1.);
	TComplex vN = FFactor::W(k0,t0v,a[1].val,1.);

	/*// Mass and width of Rho4P are parameters
	mV[4] = a[14].val;
	wV[4] = a[15].val;
	mV2[4] = (mV[4]-kI*wV[4]/2.)*(mV[4]-kI*wV[4]/2.);*/	
	
	double sign, mt;
	for (int i = 0; i < FF[1].mesons; i++) {
		mt = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if ( mt < a[1].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = FFactor::W(mwV2[i],t0v,a[1].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,4);	
	suma = normA*(FF[1].nor*mul[1]*mul[2] +
	  (mul[0]*mul[2]*(sub[2]-sub[0])/(sub[2]-sub[1]) +
	   mul[1]*mul[0]*(sub[1]-sub[0])/(sub[1]-sub[2]) -
	   mul[1]*mul[2])*a[8].val);
	return suma;
}

/** 0 = Rho, 1 = Rho1P, 2 = Rho2P, 3 = Rho3P */

TComplex FFactor::VectorTwo (TComplex t)
{
	TComplex v = FFactor::W(t,t0v,a[3].val,1.);
	TComplex vN = FFactor::W(k0,t0v,a[3].val,1.);

	/*// Mass and width of Rho4P are parameters
	mV[4] = a[14].val;
	wV[4] = a[15].val;	
	mV2[4] = (mV[4]-kI*wV[4]/2.)*(mV[4]-kI*wV[4]/2.);*/	
	
	double sign, mt;
	for (int i = 0; i < FF[3].mesons; i++) {
		mt = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if ( mt < a[3].val ) {
			sign = 1.;
		}
		else {
			sign = -1.;
		}
		vM[i] = FFactor::W(mwV2[i],t0v,a[3].val,sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,6);
	suma = normA*FF[3].nor*mul[0]*mul[1]*mul[2];
	return suma;
}

TComplex FFactor::GE1 (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a+b+t/(4.*massH1*massH1)*(c+d);
}

TComplex FFactor::GM1 (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a+b+c+d;
}

TComplex FFactor::GE2 (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a-b+t/(4.*massH2*massH2)*(c-d);
}

TComplex FFactor::GM2 (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a-b+c-d;
}

double FFactor::GE1N ( void ) { return FF[0].nor + FF[1].nor; }
double FFactor::GM1N ( void ) { return FF[0].nor + FF[1].nor + FF[2].nor + FF[3].nor; }
double FFactor::GE2N ( void ) { return FF[0].nor - FF[1].nor; }
double FFactor::GM2N ( void ) { return FF[0].nor - FF[1].nor + FF[2].nor - FF[3].nor; }

double FFactor::AbsGE1 (TComplex t)
{
	return FFactor::GE1(t).Rho();
}

double FFactor::AbsGM1 (TComplex t)
{
	return FFactor::GM1(t).Rho();
}

double FFactor::AbsGE2 (TComplex t)
{
	return FFactor::GE2(t).Rho();
}

double FFactor::AbsGM2 (TComplex t)
{
	return FFactor::GM2(t).Rho();
}

double FFactor::SigmaTotal1 ( const TComplex &t )
{
	TComplex a = a.Abs(FFactor::GM1(t));
	TComplex b = b.Abs(FFactor::GE1(t));
	TComplex k = 4./3.*TMath::Pi()*alpha*alpha/t;
	TComplex r = r.Sqrt(1.-4*massH1*massH1/t);
	TComplex v = k*r*(a*a+2.*massH1*massH1/t*b*b);
	return v;
}

double FFactor::SigmaTotal2 ( const TComplex &t )
{
	TComplex a = a.Abs(FFactor::GM2(t));
	TComplex b = b.Abs(FFactor::GE2(t));
	TComplex k = 4./3.*TMath::Pi()*alpha*alpha/t;
	TComplex r = r.Sqrt(1.-4*massH2*massH2/t);
	TComplex v = k*r*(a*a+2.*massH2*massH2/t*b*b);
	return v;
}

TComplex FFactor::TypeDefVal ( const int &type, const TComplex &t, const double &theta, const double &energy )
{
	TComplex val,nE,nM;
	double tau,eps,th,Es,Mott,DipFF,InvDipFF;

	/** "protonElectric" */
	if ((type == 1) || (type == 2)) { val = FFactor::AbsGE1(t);	}
	/** "protonMagnetic" */
	if ((type == 3) || (type == 4)) { val = FFactor::AbsGM1(t);	}
	/** "neutronElectric" */
	if ((type == 5) || (type == 6)) { val = FFactor::AbsGE2(t);	}
	/** "neutronMagnetic" */
	if (type == 7) { val = FFactor::GM2(t); }
	if (type == 8) { val = FFactor::AbsGM2(t); }
	/** "protonRatios" */
	if (type == 9) {
		nE = FFactor::AbsGE1(t);
		nM = FFactor::AbsGM1(t);
		val = val.Abs((1.+ammP)*nE/nM);
	}
	/** "neutronRatios" */
	if (type == 10) {
		nE = FFactor::AbsGE2(t);
		nM = FFactor::AbsGM2(t);
		val = val.Abs(ammN*nE/nM);
	}
	/** "Mainz Ratios" */
	if (type >= 100) {		
		th = theta/180.*TMath::Pi();
		nE = FFactor::AbsGE1(t);
		nM = FFactor::AbsGM1(t);
		tau = t/(4*massP*massP);
		eps = 1./(1.+2.*(1.+tau)*tan(th/2.)*tan(th/2.));
		Es = energy/(1.+energy/massP*(1.-cos(th)));
		InvDipFF = (1.-t/0.71)*(1.-t/0.71);
		DipFF = 1./InvDipFF;
		val = (eps*nE*nE+tau*nM*nM)/(eps+tau*(1.+ammP)*(1.+ammP))*InvDipFF*InvDipFF;
	}
	return val;
}

double FFactor::RadiusE1 ( const double step )
{
	if (step > t0v) {
		std::cout << "\n> RadiusE1: Warning! Step = " << step << " must be lower than t0v = " << t0v << std::endl;
	}
	TComplex b = FFactor::GE1(step);
	TComplex a = FFactor::GE1(-step);
	TComplex v2 = 3.*(b-a)/step*0.1*hTransC2;
	double FN = 1.;	
	std::cout << "> " << FN << " " << FFactor::GE1(0.) << std::endl;	
	return v2.Re();
}

double FFactor::RadiusE2 ( const double step )
{
	TComplex b = FFactor::GE2(step);
	TComplex a = FFactor::GE2(-step);
	TComplex v2 = 3.*(b-a)/step*0.1*hTransC2;

	double FN = 1.;
	if ( (particleType == 3) || (particleType == 5) ) { FN = -1.; } // norm. factor
	std::cout << "> " << FN << " " << FFactor::GE2(0.) << std::endl;
	return v2.Re()/FN;
}

double FFactor::RadiusM1 ( const double step )
{
	TComplex b = FFactor::GM1(step);
	TComplex a = FFactor::GM1(-step);
	TComplex v2 = 3.*(b-a)/step*0.1*hTransC2;
	double FN = FFactor::GM1N();
	std::cout << "> " << b.Re() << std::endl;
	std::cout << "> " << a.Re() << std::endl;
	std::cout << "> " << b-a << std::endl;
	std::cout << "> " << FN << " " << FFactor::GM1(0.) << std::endl;
	return v2.Re()/FN;
}

double FFactor::RadiusM2 ( const double step )
{
	TComplex b = FFactor::GM2(step);
	TComplex a = FFactor::GM2(-step);
	TComplex v2 = 3.*(b-a)/step*0.1*hTransC2;
	double FN = FFactor::GM2N();
	std::cout << "> " << b.Re() << std::endl;
	std::cout << "> " << a.Re() << std::endl;
	std::cout << "> " << b-a << std::endl;
	std::cout << "> " << FN << " " << FFactor::GM2(0.) << std::endl;
	return v2.Re()/FN;
}

// Options: all
void FFactor::CheckFormFactor ( const char* nameString, double sD )
{
	std::cout << "\n>> Standard checks for `" << nameString << "' form factors is empty!" << std::endl;
}

void FFactor::CheckParameters ()
{
	std::cout << "\n>> Standard checks for parameters:" << std::endl;

	handSome = false;
	double mH = 0;
	double mL = 0;
	int j;
	for (int i = 0; i < 4; ++i) {
		if ( (v[i].val < v[i].down) || (v[i].val > v[i].up) ) {
			std::cout << "> CheckParameters: Error!" << std::endl;
			std::cout << v[i].name << " ";
			std::cout.width(12);
			std::cout << v[i].down << " ";
			std::cout.width(12);
			std::cout << v[i].val << " ";
			std::cout.width(12);
			std::cout << v[i].up << std::endl;
			handSome = true;
		}		
	}

	if (handSome == false) {
		std::cout << "[ OK ] ... parameters" << std::endl;
	} else {
		std::cout << "> Change value or limits for parameter(s)" << std::endl;
		std::cout << "S: ";
		for (int j = 0; j < FF[0].mesons; ++j) {
			std::cout << mS[j]*mS[j]-wS[j]*wS[j]/4. << " ";
		}
		std::cout << std::endl;

		std::cout << "V: ";
		for (int j = 0; j < FF[1].mesons; ++j) {
			std::cout << mV[j]*mV[j]-wV[j]*wV[j]/4. << " ";
		}
		std::cout << std::endl;
	}

}

void FFactor::CheckCovMatrix ()
{
	std::cout << "\n>> Standard check for covariance matrix:" << std::endl;

	handSome = false;
	double h = 0.0001;
	double df,aa,bb;
	for (int i = 0; i < modelPar; ++i) {
		aa = sqrt(cov(i,i));
		bb = v[i].err;
		df = fabs(aa-bb);
		if ( df >= h ) {
			std::cout << "\n> CheckCovMatrix: Error! Parameter close to limits:     " << v[i].name << std::endl;
			std::cout << "> Error value : Sqrt of covariance element : Difference" << std::endl;
			std::cout.width(12);
			std::cout << bb << "   ";
			std::cout.width(12);
			std::cout << aa << "   ";
			std::cout.width(12);
			std::cout << df << std::endl;
			handSome = true;
		}		
	}

	if (handSome == false) {
		std::cout << "> [OK] ... covariance matrix" << std::endl;
	}	

}

	}  // namespace Minuit2

}  // namespace ROOT
