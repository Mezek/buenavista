/**
 * $Date: 2013-06-19 13:54:42 +0200 (Wed, 19 Jun 2013) $
 * $Revision: 379 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/modules/Nucleon3GH.cpp $
 * $Id: Nucleon3GH.cpp 379 2013-06-19 11:54:42Z bartos $
 *
 * @file
 * @brief	U&A form factor for nucleons &mdash; 3 vector meson resonances, v. hyperon.
 */

#include "NucleonUamH.h"

/**
 * Namespace for some useful functions.
 */

namespace shortCut {

/// Compare difference of two numbers with a small delta.

int howDiff ( double a, double b, double stdDiff )
{
	double p = fabs(a-b);
	//std::cout << ">> " << p << std::endl;
	if (p <= stdDiff) { return 0; } else { return 1; };
}

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

FFactor::FFactor (): a(globalFF), v(globalFF), particleType(1)
{
	FFactor::SetParticle(particleType);
}

FFactor::FFactor ( int myParticle ): a(globalFF), v(globalFF), particleType(myParticle)
{
	FFactor::SetParticle(particleType);
}

void FFactor::SetParticle ( int myParticle )
{
	//std::cout << "\n>> FFactor has type: " << myParticle << std::endl;
	particleType = myParticle;
	
	modelPar = globalFF;
	numberOfParameters = modelPar;
	handSome = false;

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
		mS2[i] = mS[i]*mS[i];
		mV2[i] = mV[i]*mV[i];
	}	

	// izo :: name, nor, mesons, tresh

	// Names
	FF[0].name = "F1s";
	FF[1].name = "F1v";
	FF[2].name = "F2s";
	FF[3].name = "F2v";

	// Norms
	while ( HyperonsCheckRange(particleType) != 1 ) {
		particleType = HyperonsGetType();
		std::cout << ">> FFactor has type: " << particleType << std::endl;
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

	// Mesons under threshold
	FF[0].tresh = 3;
	FF[1].tresh = 2;
	FF[2].tresh = 3;
	FF[3].tresh = 2;
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
			++k;
		}
		myDataFile.close();
		if (N != modelPar) {
			std::cout << "> LoadParameters: Error!" << std::endl;
			std::cout << "> Number of parameters: " << N << " in '" << ds << "' differs from model: " << modelPar << "!" << std::endl;
		}
	}
	else std::cerr << ">> LoadParameters: Error! Unable to open parametric file: '" << ds << "'!" << std::endl;

	FFactor::FixParameters();
	// debug: if ( particleType == 1 ) { particleType = 5; }	
	if ( particleType != 1 ) { FFactor::TransformSU3(); }
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
		std::cout.precision(10);
		std::cout.width(17);
		std::cout << this->v[i].val << "   +/- ";
		std::cout.width(17);
		std::cout << this->v[i].err << std::endl;
	}

	std::cout << "\n>> Norms for the components:" << std::endl;
	for (int i = 0; i < 4; ++i) {
		std::cout << FF[i].name << ": nor[" << i << "] = " << FF[i].nor << std::endl;
	}

	// Signs of parameters
	std::cout.width(10);
	std::cout.precision(7);
	std::cout << "\n>> Signs of the parameters:" << std::endl;
	std::cout << "\t" << "Om\t" << "Phi\t" << "Om1\t" << "Ph1\t"<< "Rh\t" << std::endl;
	std::cout << "F1s\t" << shortCut::signum(this->v[4].val) << "\t" <<
		shortCut::signum(this->v[5].val) << "\t" <<
		shortCut::signum(this->v[6].val) << "\t" <<
		shortCut::signum(this->v[7].val) << std::endl;
	std::cout << "F2s\t" << shortCut::signum(this->v[9].val) << "\t" <<
		shortCut::signum(this->v[10].val) << "\t" <<
		shortCut::signum(this->v[11].val) << std::endl;
	std::cout << "F1v\t" << "\t" << "\t"<< "\t" << "\t" <<
		shortCut::signum(this->v[8].val) << std::endl;
	std::cout << "F2v\t" << "\t" << "\t"<< "\t" << "\t" <<
		shortCut::signum(this->v[12].val) << std::endl;

	// Mesons under/over threshold
	std::cout << "\n>> Mesons near threshold:" << std::endl;
	std::cout << "\t" << "Under\t" << "Over\t" << "All\t" << std::endl;

	for (int i = 0; i < 4; ++i) {
		std::cout << FF[i].name << "\t" << FF[i].tresh << "\t" << 
				     FF[i].mesons-FF[i].tresh << "\t" << FF[i].mesons << std::endl;
	}
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
	double th = 40.46*TMath::DegToRad();
	double k1 = TMath::Cos(th)/TMath::Sqrt(2.);
	double k2 = TMath::Sin(th)/TMath::Sqrt(3.);
	double k3 = TMath::Sin(th)/TMath::Sqrt(2.);
	double k4 = TMath::Cos(th)/TMath::Sqrt(3.);

	double fuO[3], fuP[3], fuR[3]; // universal coupling constants
	fuO[0] = 17.0576;
	fuP[0] = 13.4448;
	fuR[0] = 4.9569;
	fuO[1] = 47.5897;
	fuP[1] = 37.0510;
	fuR[1] = 13.6455;
	fuO[2] = 48.3651;
	fuP[2] = 0.;
	fuR[2] = 22.5275;

	double fO[2][3], fP[2][3], fR[2][3];
	fO[0][0] = a[4].val*fuO[0];
	fP[0][0] = a[5].val*fuP[0];
	fR[0][0] = a[8].val*fuR[0];
	fO[1][0] = a[9].val*fuO[0];
	fP[1][0] = a[10].val*fuP[0];
	fR[1][0] = a[12].val*fuR[0];

	fO[0][1] = a[6].val*fuO[1];
	fP[0][1] = a[7].val*fuP[1];
	fR[0][1] = (FF[1].nor*Duo(mV2[2],mV2[1]) -
		Trio(mV2[2],mV2[0],mV2[1])*a[8].val)*fuR[1];
	fO[1][1] = a[11].val*fuO[1];
	fP[1][1] = (FF[2].nor*Duo(mS2[4],mS2[3])*Duo(mS2[5],mS2[3]) -
		a[9].val*Trio(mS2[4],mS2[0],mS2[3])*Trio(mS2[5],mS2[0],mS2[3]) -
		a[10].val*Trio(mS2[4],mS2[1],mS2[3])*Trio(mS2[5],mS2[1],mS2[3]) -
		a[11].val*Trio(mS2[4],mS2[2],mS2[3])*Trio(mS2[5],mS2[2],mS2[3]))*fuP[1];
	fR[1][1] = (FF[3].nor*Duo(mV2[2],mV2[1]) -
		Trio(mV2[2],mV2[0],mV2[1])*a[12].val)*fuR[1];

	fO[0][2] = (FF[0].nor*Duo(mS2[5],mS2[4]) -
		a[4].val*Trio(mS2[5],mS2[0],mS2[4]) -
		a[5].val*Trio(mS2[5],mS2[1],mS2[4]) -
		a[6].val*Trio(mS2[5],mS2[2],mS2[4]) -
		a[7].val*Trio(mS2[5],mS2[3],mS2[4]))*fuO[2];
	fP[0][2] = (FF[0].nor*Duo(mS2[4],mS2[5]) -
		a[4].val*Trio(mS2[4],mS2[0],mS2[5]) -
		a[5].val*Trio(mS2[4],mS2[1],mS2[5]) -
		a[6].val*Trio(mS2[4],mS2[2],mS2[5]) -
		a[7].val*Trio(mS2[4],mS2[3],mS2[5]))*fuP[2];
	fR[0][2] = (FF[1].nor*Duo(mV2[1],mV2[2]) -
		Trio(mV2[1],mV2[0],mV2[2])*a[8].val)*fuR[2];
	fO[1][2] = (FF[2].nor*Duo(mS2[3],mS2[4])*Duo(mS2[5],mS2[4]) -
		a[9].val*Trio(mS2[3],mS2[0],mS2[4])*Trio(mS2[5],mS2[0],mS2[4]) -
		a[10].val*Trio(mS2[3],mS2[1],mS2[4])*Trio(mS2[5],mS2[1],mS2[4]) -
		a[11].val*Trio(mS2[3],mS2[2],mS2[4])*Trio(mS2[5],mS2[2],mS2[4]))*fuO[2];
	fP[1][2] = (FF[2].nor*Duo(mS2[3],mS2[5])*Duo(mS2[4],mS2[5]) -
		a[9].val*Trio(mS2[3],mS2[0],mS2[5])*Trio(mS2[4],mS2[0],mS2[5]) -
		a[10].val*Trio(mS2[3],mS2[1],mS2[5])*Trio(mS2[4],mS2[1],mS2[5]) -
		a[11].val*Trio(mS2[3],mS2[2],mS2[5])*Trio(mS2[4],mS2[2],mS2[5]))*fuP[2];
	fR[1][2] = (FF[3].nor*Duo(mV2[1],mV2[2]) -
		Trio(mV2[1],mV2[0],mV2[2])*a[12].val)*fuR[2];

	// Inverse equations
	double fD[2][3], fF[2][3], fS[2][3];
	double fOp[2][3], fPp[2][3], fRp[2][3];
	for (int j = 0; j < 3; ++j) {
		for (int i = 0; i < 2; ++i) {
			fD[i][j] = ((6.*k1*k4+6.*k2*k3)*fR[i][j] -
				2.*k1*fP[i][j]+2.*k3*fO[i][j])/(4.*k1*k4+k3*k3+3.*k2*k3);
			fF[i][j] = ((2.*k1*k4+2.*k3*k3)*fR[i][j] +
				2.*k1*fP[i][j]-2.*k3*fO[i][j])/(4.*k1*k4+k3*k3+3.*k2*k3);
			fS[i][j] = -((3.*k3*k4-3*k2*k4)*fR[i][j]-k3*fP[i][j] -
				3.*k2*fP[i][j]-4.*k4*fO[i][j])/(4.*k1*k4+k3*k3+3.*k2*k3);
			//std::cout << i << " " << fD[i][j] << " " << fF[i][j] << " " << fS[i][j] << std::endl;

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
	a[11].val = fOp[1][1]/fuO[1];
	a[12].val = fRp[1][0]/fuR[0];
	FFactor::FixParameters();	
	std::cout << ">> SU(3) transformation of coupling constants done." << std::endl;
	
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

/** 0 = Om, 1 = Phi, 2 = Om1P, 3 = Phi1P, 4 = Om2P, 5 = Phi2P */

TComplex FFactor::ScalarOne (TComplex t)
{
	TComplex v = FFactor::W(t,t0s,a[0].val,1.);
	TComplex vN = FFactor::W(k0,t0s,a[0].val,1.);
	//std::cout << ">> ScalarOne: "<< a[0] << " " << v << " " << vN << std::endl;
	
	double sign,threshold;
	for (int i = 0; i < FF[0].mesons; i++) {
		threshold = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		//if (threshold > a[0]) {
		if ( i < FF[0].tresh ) {
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
	suma = normA*(FF[0].nor*mul[0]*mul[1] +
	  (mul[2]*mul[1]*(sub[1]-sub[2])/(sub[1]-sub[0]) +
	   mul[2]*mul[0]*(sub[0]-sub[2])/(sub[0]-sub[1]) -
	   mul[0]*mul[1])*a[4].val +
	  (mul[3]*mul[1]*(sub[1]-sub[3])/(sub[1]-sub[0]) +
	   mul[3]*mul[0]*(sub[0]-sub[3])/(sub[0]-sub[1]) -
	   mul[0]*mul[1])*a[5].val +
	  (mul[4]*mul[1]*(sub[1]-sub[4])/(sub[1]-sub[0]) +
	   mul[4]*mul[0]*(sub[0]-sub[4])/(sub[0]-sub[1]) -
	   mul[0]*mul[1])*a[6].val +
	  (mul[5]*mul[1]*(sub[1]-sub[5])/(sub[1]-sub[0]) +
	   mul[5]*mul[0]*(sub[0]-sub[5])/(sub[0]-sub[1]) -
	   mul[0]*mul[1])*a[7].val);

	//std::cout << ">> ScalarOne: "<< norm << " " << FF[0].nor << std::endl;
	return suma;
}

/** 0 = Om, 1 = Phi, 2 = Om1P, 3 = Phi1P, 4 = Om2P, 5 = Phi2P */

TComplex FFactor::ScalarTwo (TComplex t)
{
	TComplex v = FFactor::W(t,t0s,a[2].val,1.);
	TComplex vN = FFactor::W(k0,t0s,a[2].val,1.);
	
	double sign,threshold;
	for (int i = 0; i < FF[2].mesons; i++) {
		threshold = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		//if (threshold > a[2]) {
		if ( i < FF[2].tresh ) {
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
	suma = normA*(FF[2].nor*mul[0]*mul[1]*mul[2] +
	  (mul[3]*mul[1]*mul[2]*(sub[1]-sub[3])/(sub[1]-sub[0])*(sub[2]-sub[3])/(sub[2]-sub[0]) +
	   mul[3]*mul[0]*mul[2]*(sub[0]-sub[3])/(sub[0]-sub[1])*(sub[2]-sub[3])/(sub[2]-sub[1]) +
	   mul[3]*mul[0]*mul[1]*(sub[0]-sub[3])/(sub[0]-sub[2])*(sub[1]-sub[3])/(sub[1]-sub[2]) -
	   mul[0]*mul[1]*mul[2])*a[9].val +
	  (mul[4]*mul[1]*mul[2]*(sub[1]-sub[4])/(sub[1]-sub[0])*(sub[2]-sub[4])/(sub[2]-sub[0]) +
	   mul[4]*mul[0]*mul[2]*(sub[0]-sub[4])/(sub[0]-sub[1])*(sub[2]-sub[4])/(sub[2]-sub[1]) +
	   mul[4]*mul[0]*mul[1]*(sub[0]-sub[4])/(sub[0]-sub[2])*(sub[1]-sub[4])/(sub[1]-sub[2]) -
	   mul[0]*mul[1]*mul[2])*a[10].val +
	  (mul[5]*mul[1]*mul[2]*(sub[1]-sub[5])/(sub[1]-sub[0])*(sub[2]-sub[5])/(sub[2]-sub[0]) +
	   mul[5]*mul[0]*mul[2]*(sub[0]-sub[5])/(sub[0]-sub[1])*(sub[2]-sub[5])/(sub[2]-sub[1]) +
	   mul[5]*mul[0]*mul[1]*(sub[0]-sub[5])/(sub[0]-sub[2])*(sub[1]-sub[5])/(sub[1]-sub[2]) -
	   mul[0]*mul[1]*mul[2])*a[11].val);

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
	
	double sign,threshold;
	for (int i = 0; i < FF[1].mesons; i++) {
		threshold = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		//if (threshold > a[1]) {
		if ( i < FF[1].tresh ) {
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
	suma = normA*(FF[1].nor*mul[0]*mul[1] +
	  (mul[2]*mul[1]*(sub[1]-sub[2])/(sub[1]-sub[0]) +
	   mul[0]*mul[2]*(sub[0]-sub[2])/(sub[0]-sub[1]) -
	   mul[0]*mul[1])*a[8].val);
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
	
	double sign,threshold;
	for (int i = 0; i < FF[3].mesons; i++) {
		threshold = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		//if (threshold > a[3]) {
		if ( i < FF[3].tresh ) {
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
	suma = normA*(FF[3].nor*mul[0]*mul[1] +
	  (mul[2]*mul[1]*(sub[1]-sub[2])/(sub[1]-sub[0]) +
	   mul[0]*mul[2]*(sub[0]-sub[2])/(sub[0]-sub[1]) -
	   mul[0]*mul[1])*a[12].val);
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
		if ((i == 0) || (i == 2)) {
			j = FF[i].tresh;
			mH = mS[j]*mS[j]-wS[j]*wS[j]/4.;
			mL = mS[j-1]*mS[j-1]-wS[j-1]*wS[j-1]/4.;
		} else {
			mH = mV[j]*mV[j]-wV[j]*wV[j]/4.;
			mH = mV[j-1]*mV[j-1]-wV[j-1]*wV[j-1]/4.;
		}
		if ( v[i].val > mH ) {
			std::cout << "> CheckParameters: Error!" << std::endl;
			std::cout << "> Change value or limits for parameter: " << v[i].name << std::endl;
			std::cout << "> Value:            " << v[i].val << std::endl;
			std::cout << "> Higher limit:     " << mH << std::endl;
			handSome = true;
		}		
		if ( v[i].val < mL ) {
			std::cout << "> CheckParameters: Error!" << std::endl;
			std::cout << "> Change value or limits for parameter: " << v[i].name << std::endl;
			std::cout << "> Value:            " << v[i].val << std::endl;
			std::cout << "> Lower limit:      " << mL << std::endl;
			handSome = true;
		}		
	}

	if (handSome == false) {
		std::cout << "> [OK] ... parameters" << std::endl;
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
