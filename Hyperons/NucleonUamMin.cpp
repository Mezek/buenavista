/**
 * $Date: 2013-05-22 10:09:53 +0200 (Wed, 22 May 2013) $
 * $Revision: 349 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/NucleonUamMin.cpp $
 * $Id: NucleonUamMin.cpp 349 2013-05-22 08:09:53Z bartos $
 *
 * @file
 * @brief	UA form factor for nucleons - minimal model - using SU(3)
 */

#include <NucleonUamH.h>

namespace ROOT {

	namespace Minuit2 {
		
FFactor::FFactor ( std::size_t size, const int& myParticle ): a(size), val(size), particleType(myParticle)
{ 
	t(1.,0.00001);
	
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
	
	for (int i=0; i<6; i++) {
		mS2[i] = (mS[i]-kI*wS[i]/2.)*(mS[i]-kI*wS[i]/2.);
		mV2[i] = (mV[i]-kI*wV[i]/2.)*(mV[i]-kI*wV[i]/2.);
	}

	modelPar = a.size();
	numberOfParameters = modelPar;

	while ( HyperonsCheckRange(particleType) != 1 ) {
		particleType = HyperonsGetType();
	}

	switch (particleType) {
		case 1:
			nor[0] = 0.5;
			nor[1] = 0.5*(ammP+ammN);
			nor[2] = 0.5;
			nor[3] = 0.5*(ammP-ammN);
			massH1 = massP;
			massH2 = massN;
			break;
		case 2:
			nor[0] = 0.;
			nor[1] = muL;
			nor[2] = 0.;
			nor[3] = 0.;
			massH1 = massL;
			massH2 = 0.;
			break;
		case 3:
			nor[0] = 0.;
			nor[1] = 0.5*(muSp+muSm);
			nor[2] = 1.;
			nor[3] = 0.5*(muSp-muSm);
			massH1 = massSp;
			massH2 = massSm;
			break;
		case 4:
			nor[0] = 0.;
			nor[1] = 0.5*(muSp+muSm);
			nor[2] = 1.;
			nor[3] = 0.5*(muSp-muSm);
			massH1 = massS0;
			massH2 = 0.;
			break;			
		case 5:
			nor[0] = -0.5;
			nor[1] = 0.5*(muX0+muXm);
			nor[2] = 0.5;
			nor[3] = 0.5*(muX0-muXm);
			massH1 = massX0;
			massH2 = massXm;			
			break;
		default:
			std::cout << ">> Error: No particle type defined!" << std::endl;
	}
}

void FFactor::LoadParameters (char* ds)
{
	ifstream myDataFile (ds);
	int num = 0;
	if (myDataFile.is_open()) {
		while (myDataFile.peek() != EOF) {
			myDataFile >> a[num];
			this->val[num]=a[num];
			//std::cout << num << " " << a[num] << " " << val[num] << std::endl;
			num++;
		}
		myDataFile.close();
		if (num-1 != modelPar) {
			std::cout << ">> Warning: Number of parameters: " << num-1 << " in '" << ds << "' differs from model: " << modelPar << "!" << std::endl;
		}
	}
	else std::cerr << ">> Error: Unable to open parametric file: '" << ds << "'!" << std::endl;
}

void FFactor::SetParameters (const std::vector<double>& par)
{
	for (int i=0; i<modelPar; i++) 
	{
		a[i] = par[i];
		val[i] = par[i];		
		//std::cout << i << ". " << a[i] << " " << val[i] << std::endl;
	}
}

void FFactor::FixParameters () {}

void FFactor::PrintParameters ()
{
	std::cout << ">> Actual pararameters of FFactor:" << std::endl;
	for (int i=0; i<modelPar; i++) 
	{
		//std::cout << i+1 << ". " << a[i] << std::endl;
		std::cout << i+1 << ".: " << this->val[i] << std::endl;
	}

	for (int i=0; i<4; i++) {
		std::cout << "nor[" << i << "] = " << nor[i] << std::endl;
	}
}

void FFactor::TransformSU3 ()
{
	// Calculate SU(3) relations
	double th = 39.83*TMath::DegToRad();
	double k1 = TMath::Cos(th)/TMath::Sqrt(2.);
	double k2 = TMath::Sin(th)/TMath::Sqrt(3.);
	double k3 = TMath::Sin(th)/TMath::Sqrt(2.);
	double k4 = TMath::Cos(th)/TMath::Sqrt(3.);

	double fO[2], fP[2], fR[2];
	double fD[2], fF[2], fS[2];

	double fdO = 17.06;
	double fdP = 1.27;
	double fdR = 4.96;

	fO[0] = val[4]*fdO;
	fP[0] = val[5]*fdP;
	fR[0] = val[6]*fdR;
	fO[1] = val[7]*fdO;
	fP[1] = val[8]*fdP;
	fR[1] = val[9]*fdR;

	// Inverse equations
	for (int i = 0; i < 2; i++) {
		fD[i] = ((6.*k1*k4+6.*k2*k3)*fR[i]-2.*k1*fP[i]+2.*k3*fO[i])/(4.*k1*k4+k3*k3+3.*k2*k3);
		fF[i] = ((2.*k1*k4+2.*k3*k3)*fR[i]+2.*k1*fP[i]-2.*k3*fO[i])/(4.*k1*k4+k3*k3+3.*k2*k3);
		fS[i] = -((3.*k3*k4-3*k2*k4)*fR[i]-k3*fP[i]-3.*k2*fP[i]-4.*k4*fO[i])/(4.*k1*k4+k3*k3+3.*k2*k3);
		//std::cout << i << " " << fD[i] << " " << fF[i] << " " << fS[i] << std::endl;
	}

	double fOlam[2], fOsig[2], fOksi[2];
	double fPlam[2], fPsig[2], fPksi[2];
	double fRsig[2], fRksi[2];

	for (int i = 0; i < 2; i++) {
		fOlam[i] = k1*fS[i] + k2*fD[i];
		fPlam[i] = k3*fS[i] - k4*fD[i];

		fOsig[i] = k1*fS[i] - k3*fD[i];
		fPsig[i] = k3*fS[i] + k4*fD[i];
		fRsig[i] = -fF[i];

		fOksi[i] = k1*fS[i] + 3./2.*k2*fF[i] + k2/2.*fD[i];
		fPksi[i] = k3*fS[i] - 3./2.*k4*fF[i] - k4/2.*fD[i];
		fRksi[i] = (fD[i] - fF[i])/2.;
	}

	val[4] = fOsig[0]/fdO;
	val[5] = fPsig[0]/fdP; 
	val[6] = fRsig[0]/fdR;
	val[7] = fOsig[1]/fdO;
	val[8] = fPsig[1]/fdP;
	val[9] = fRsig[1]/fdR;
}

// sign = +1. : below threshold - Low
// sign = -1. : over threshold - High
// It coressponds to change V(t) to -1./V(t)

TComplex FFactor::W (const TComplex &t, const TComplex &t0, const TComplex &tin, const double sign)
{
	TComplex cQ,cQin,cK,cZ,cRes;

	cQ = cQ.Sqrt((t-t0)/t0);
	cQin = cQin.Sqrt((tin-t0)/t0);

	if (cQ.Re() > tin.Re()) {
		cQ = cQ + 0.000001*kI;
		//std::cout << "cQ: " << a[0] << " " << t << " " << cQ << std::endl;
	}

	cK = cK.Sqrt(cQin+cQ);
	cZ = cZ.Sqrt(cQin-cQ);
	cRes = kI*(cK-sign*cZ)/(cK+sign*cZ);
	return cRes;
}

TComplex eL (const TComplex &a, const TComplex &b, const TComplex &c, const TComplex &sc)
{
	TComplex f = (b-c)/(a-c)*(b-sc)/(a-sc)*(b-1./c)/(a-1./c)*(b-1./sc)/(a-1./sc);
	return f;
}

TComplex eH (const TComplex &a, const TComplex &b, const TComplex &c, const TComplex &sc)
{
	TComplex f = (b-c)/(a-c)*(b-sc)/(a-sc)*(b+c)/(a+c)*(b+sc)/(a+sc);
	return f;
}

TComplex sIL (const TComplex &b, const TComplex &c, const TComplex &sc)
{
	TComplex f = -1.*(b-c)*(b-sc)/(c-1./c)*(b-1./c)*(b-1./sc)/(sc-1./sc);
	return f;
}

TComplex sIH (const TComplex &b, const TComplex &c, const TComplex &sc)
{
	TComplex f = -1.*(b-c)*(b-sc)/(c-1./c)*(b+c)*(b+sc)/(sc-1./sc);
	return f;
}

// 0 = Om, 1 = Phi, 2 = Om1P, 3 = Phi1P, 4 = Om2P, 5 = Phi2P

TComplex FFactor::ScalarOne (TComplex t)
{
	TComplex v = FFactor::W(t,t0s,a[0],1.);
	TComplex vN = FFactor::W(k0,t0s,a[0],1.);
	//std::cout << ">> ScalarE: "<< a[0] << " " << v << " " << vN << std::endl;
	
	double threshold;
	
	for (int i=0; i<5; i++) {
		threshold = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if (threshold > a[0]) {
			vM[i] = FFactor::W(mS2[i],t0s,a[0],-1.);
			vMc[i] = vMc[i].Conjugate(vM[i]);
			mul[i] = eH(v,vN,vM[i],vMc[i]);
			//sub[i] = sIH(vN,vM[i],vMc[i]);
		}
		else {
			vM[i] = FFactor::W(mS2[i],t0s,a[0],1.);
			vMc[i] = vMc[i].Conjugate(vM[i]);
			mul[i] = eL(v,vN,vM[i],vMc[i]);
			//sub[i] = sIL(vN,vM[i],vMc[i]);
		}
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,4);
	suma = normA*((mul[0]-mul[2])*a[4] + (mul[1]-mul[2])*a[5] + mul[2]*nor[0]);
	return suma;
}

TComplex FFactor::ScalarTwo (TComplex t)
{
	TComplex v = FFactor::W(t,t0s,a[2],1.);
	TComplex vN = FFactor::W(k0,t0s,a[2],1.);
	
	double threshold;
	
	for (int i=0; i<5; i++) {
		threshold = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if (threshold > a[2]) {
			vM[i] = FFactor::W(mS2[i],t0s,a[2],-1.);
			vMc[i] = vMc[i].Conjugate(vM[i]);
			mul[i] = eH(v,vN,vM[i],vMc[i]);
			//sub[i] = sIH(vN,vM[i],vMc[i]);
		}
		else {
			vM[i] = FFactor::W(mS2[i],t0s,a[2],1.);
			vMc[i] = vMc[i].Conjugate(vM[i]);
			mul[i] = eL(v,vN,vM[i],vMc[i]);
			//sub[i] = sIL(vN,vM[i],vMc[i]);
		}
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,6);
	suma = normA*((mul[0]-mul[2])*a[7] + (mul[1]-mul[2])*a[8] + mul[2]*nor[1]);
	return suma;
}

// 0 = Rho, 1 = Rho1P, 2 = Rho2P, 3 = Rho3P

TComplex FFactor::VectorOne (TComplex t)
{
	TComplex v = FFactor::W(t,t0v,a[1],1.);
	TComplex vN = FFactor::W(k0,t0v,a[1],1.);
	
	double threshold;
	for (int i=0; i<3; i++) {
		threshold = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if (threshold > a[1]) {
			vM[i] = FFactor::W(mV2[i],t0v,a[1],-1.);
			vMc[i] = vMc[i].Conjugate(vM[i]);
			mul[i] = eH(v,vN,vM[i],vMc[i]);
			//sub[i] = sIH(vN,vM[i],vMc[i]);
		}
		else {
			vM[i] = FFactor::W(mV2[i],t0v,a[1],1.);
			vMc[i] = vMc[i].Conjugate(vM[i]);
			mul[i] = eL(v,vN,vM[i],vMc[i]);
			//sub[i] = sIL(vN,vM[i],vMc[i]);
		}
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,4);
	suma = normA*((mul[0]-mul[1])*a[6] + mul[1]*nor[2]);
	return suma;
}

TComplex FFactor::VectorTwo (TComplex t)
{
	TComplex v = FFactor::W(t,t0v,a[3],1.);
	TComplex vN = FFactor::W(k0,t0v,a[3],1.);
	
	double sign,threshold;
	for (int i=0; i<4; i++) {
		threshold = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if (threshold > a[3]) {
			vM[i] = FFactor::W(mV2[i],t0v,a[3],-1.);
			vMc[i] = vMc[i].Conjugate(vM[i]);
			mul[i] = eH(v,vN,vM[i],vMc[i]);
			//sub[i] = sIH(vN,vM[i],vMc[i]);
		}
		else {
			vM[i] = FFactor::W(mV2[i],t0v,a[3],1.);
			vMc[i] = vMc[i].Conjugate(vM[i]);
			mul[i] = eL(v,vN,vM[i],vMc[i]);
			//sub[i] = sIL(vN,vM[i],vMc[i]);
		}
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,6);
	suma = normA*((mul[0]-mul[1])*a[9] + mul[1]*nor[3]);
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

double FFactor::SigmaTotal1 ( const TComplex &t )
{
	TComplex a = a.Abs(FFactor::GM1(t));
	TComplex b = b.Abs(FFactor::GE1(t));
	TComplex k = 4./3.*TMath::Pi()*alpha*alpha/t;
	TComplex r = r.Sqrt(1.-4*massH1*massH1/t);
	TComplex v =  k*r*(a*a+2.*massH1*massH1/t*b*b);
	return v;
}

double FFactor::SigmaTotal2 ( const TComplex &t )
{
	TComplex a = a.Abs(FFactor::GM2(t));
	TComplex b = b.Abs(FFactor::GE2(t));
	TComplex k = 4./3.*TMath::Pi()*alpha*alpha/t;
	TComplex r = r.Sqrt(1.-4*massH2*massH2/t);
	TComplex v =  k*r*(a*a+2.*massH2*massH2/t*b*b);
	return v;
}
	}  // namespace Minuit2

}  // namespace ROOT
