/**
 * $Date: 2013-05-27 14:56:24 +0200 (Mon, 27 May 2013) $
 * $Revision: 356 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/Nucleon3G.cpp $
 * $Id: Nucleon3G.cpp 356 2013-05-27 12:56:24Z bartos $
 *
 * @file
 * @brief	U&A form factor for nucleons &mdash; 3 vector meson resonances.
 */

#include "NucleonUam.h"

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
		
FFactor::FFactor ( std::size_t size ): a(size), val(size)
{ 
	t(1.,0.00001);
	//std::cout << ">> FFactor is empty! Value of t: " << t << std::endl;
	
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
		mwS2[i] = (mS[i]-kI*wS[i]/2.)*(mS[i]-kI*wS[i]/2.);
		mwV2[i] = (mV[i]-kI*wV[i]/2.)*(mV[i]-kI*wV[i]/2.);
	}

	nor[0] = 0.5;
	nor[1] = 0.5*(ammP+ammN);
	nor[2] = 0.5;
	nor[3] = 0.5*(ammP-ammN);	

	modelPar = a.size();
	numberOfParameters = modelPar;
	handSome = false;
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
			std::cout << ">> Warning! Number of parameters: " << num-1 << " in '" << ds << "' differs from model: " << modelPar << "!" << std::endl;
		}
	}
	else std::cerr << ">> Error! Unable to open parametric file: '" << ds << "'!" << std::endl;

}

void FFactor::SetParameters (const std::vector<double>& par)
{
	for (int i=0; i<modelPar; i++) 
	{
		a[i] = par[i];
		val[i] = par[i];	
	}
}

void FFactor::PrintParameters ()
{
	std::cout << "\n>> Actual pararameters of FFactor:" << std::endl;
	for (int i=0; i<modelPar; i++) 
	{
		//std::cout << i+1 << ". " << a[i] << std::endl;
		std::cout << i+1 << ". " << this->val[i] << std::endl;
	}

	for (int i=0; i<4; i++) {
		std::cout << "nor[" << i << "] = " << nor[i] << std::endl;
	}

	// Signs of parameters

	std::cout << "\n\t" << "Om\t" << "Phi\t" << "Om1\t" << "Ph1\t"<< "Rh\t" <<std::endl;
	std::cout << "F1s\t" << shortCut::signum(this->val[4]) << "\t" <<
		shortCut::signum(this->val[5]) << "\t" <<
		shortCut::signum(this->val[6]) << "\t" <<
		shortCut::signum(this->val[7]) << std::endl;
	std::cout << "F2s\t" << shortCut::signum(this->val[9]) << "\t" <<
		shortCut::signum(this->val[10]) << "\t" <<
		shortCut::signum(this->val[11]) << std::endl;
	std::cout << "F1v\t" << "\t" << "\t"<< "\t" << "\t" <<
		shortCut::signum(this->val[8]) << std::endl;
	std::cout << "F2v\t" << "\t" << "\t"<< "\t" << "\t" << "\t" << std::endl;	
}

/*
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
	TComplex v = FFactor::W(t,t0s,a[0],1.);
	TComplex vN = FFactor::W(k0,t0s,a[0],1.);
	//std::cout << ">> ScalarE: "<< a[0] << " " << v << " " << vN << std::endl;
	
	double sign,threshold;
	for (int i=0; i<6; i++) {
		threshold = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if (threshold > a[0]) {
			sign = -1.;
		}
		else {
			sign = 1.;
		}
		vM[i] = FFactor::W(mwS2[i],t0s,a[0],sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,4);
	suma = normA*(nor[0]*mul[4]*mul[5] +
	  (mul[0]*mul[5]*(sub[5]-sub[0])/(sub[5]-sub[4]) +
	   mul[0]*mul[4]*(sub[4]-sub[0])/(sub[4]-sub[5]) -
	   mul[4]*mul[5])*a[4] +
	  (mul[1]*mul[5]*(sub[5]-sub[1])/(sub[5]-sub[4]) +
	   mul[1]*mul[4]*(sub[4]-sub[1])/(sub[4]-sub[5]) -
	   mul[4]*mul[5])*a[5] +
	  (mul[2]*mul[5]*(sub[5]-sub[2])/(sub[5]-sub[4]) +
	   mul[2]*mul[4]*(sub[4]-sub[2])/(sub[4]-sub[5]) -
	   mul[4]*mul[5])*a[6] +
	  (mul[3]*mul[5]*(sub[5]-sub[3])/(sub[5]-sub[4]) +
	   mul[3]*mul[4]*(sub[4]-sub[3])/(sub[4]-sub[5]) -
	   mul[4]*mul[5])*a[7]);

	//std::cout << ">> ScalarE: "<< norm << " " << suma << std::endl;
	return suma;
}

/** 0 = Om, 1 = Phi, 2 = Om1P, 3 = Phi1P, 4 = Om2P, 5 = Phi2P */

TComplex FFactor::ScalarTwo (TComplex t)
{
	TComplex v = FFactor::W(t,t0s,a[2],1.);
	TComplex vN = FFactor::W(k0,t0s,a[2],1.);
	
	double sign,threshold;
	for (int i=0; i<6; i++) {
		threshold = mS[i]*mS[i]-wS[i]*wS[i]/4.;
		if (threshold > a[2]) {
			sign = -1.;
		}
		else {
			sign = 1.;
		}
		vM[i] = FFactor::W(mwS2[i],t0s,a[2],sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,6);
	suma = normA*(nor[1]*mul[3]*mul[4]*mul[5] +
	  (mul[0]*mul[4]*mul[5]*(sub[4]-sub[0])/(sub[4]-sub[3])*(sub[5]-sub[0])/(sub[5]-sub[3]) +
	   mul[0]*mul[3]*mul[5]*(sub[3]-sub[0])/(sub[3]-sub[4])*(sub[5]-sub[0])/(sub[5]-sub[4]) +
	   mul[0]*mul[3]*mul[4]*(sub[3]-sub[0])/(sub[3]-sub[5])*(sub[4]-sub[0])/(sub[4]-sub[5]) -
	   mul[3]*mul[4]*mul[5])*a[9] +
	  (mul[1]*mul[4]*mul[5]*(sub[4]-sub[1])/(sub[4]-sub[3])*(sub[5]-sub[1])/(sub[5]-sub[3]) +
	   mul[1]*mul[3]*mul[5]*(sub[3]-sub[1])/(sub[3]-sub[4])*(sub[5]-sub[1])/(sub[5]-sub[4]) +
	   mul[1]*mul[3]*mul[4]*(sub[3]-sub[1])/(sub[3]-sub[5])*(sub[4]-sub[1])/(sub[4]-sub[5]) -
	   mul[3]*mul[4]*mul[5])*a[10] +
	  (mul[2]*mul[4]*mul[5]*(sub[4]-sub[2])/(sub[4]-sub[3])*(sub[5]-sub[2])/(sub[5]-sub[3]) +
	   mul[2]*mul[3]*mul[5]*(sub[3]-sub[2])/(sub[3]-sub[4])*(sub[5]-sub[2])/(sub[5]-sub[4]) +
	   mul[2]*mul[3]*mul[4]*(sub[3]-sub[2])/(sub[3]-sub[5])*(sub[4]-sub[2])/(sub[4]-sub[5]) -
	   mul[3]*mul[4]*mul[5])*a[11]);
	return suma;
}

/** 0 = Rho, 1 = Rho1P, 2 = Rho2P, 3 = Rho3P */

TComplex FFactor::VectorOne (TComplex t)
{
	TComplex v = FFactor::W(t,t0v,a[1],1.);
	TComplex vN = FFactor::W(k0,t0v,a[1],1.);

	// Mass and width of Rho4P are parameters
	mV[4] = a[14];
	wV[4] = a[15];
	mV2[4] = (mV[4]-kI*wV[4]/2.)*(mV[4]-kI*wV[4]/2.);	
	
	double sign,threshold;
	for (int i=0; i<5; i++) {
		threshold = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if (threshold > a[1]) {
			sign = -1.;
		}
		else {
			sign = 1.;
		}
		vM[i] = FFactor::W(mwV2[i],t0v,a[1],sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,4);	
	suma = normA*(nor[2]*mul[1]*mul[2] +
	  (mul[0]*mul[2]*(sub[2]-sub[0])/(sub[2]-sub[1]) +
	   mul[1]*mul[0]*(sub[1]-sub[0])/(sub[1]-sub[2]) -
	   mul[1]*mul[2])*a[8]);
	return suma;
}

/** 0 = Rho, 1 = Rho1P, 2 = Rho2P, 3 = Rho3P */

TComplex FFactor::VectorTwo (TComplex t)
{
	TComplex v = FFactor::W(t,t0v,a[3],1.);
	TComplex vN = FFactor::W(k0,t0v,a[3],1.);

	// Mass and width of Rho4P are parameters
	mV[4] = a[14];
	wV[4] = a[15];	
	mV2[4] = (mV[4]-kI*wV[4]/2.)*(mV[4]-kI*wV[4]/2.);	
	
	double sign,threshold;
	for (int i=0; i<5; i++) {
		threshold = mV[i]*mV[i]-wV[i]*wV[i]/4.;
		if (threshold > a[3]) {
			sign = -1.;
		}
		else {
			sign = 1.;
		}
		vM[i] = FFactor::W(mwV2[i],t0v,a[3],sign);
		vMc[i] = vMc[i].Conjugate(vM[i]);
		mul[i] = eL(v,vN,vM[i],vMc[i],sign);
		sub[i] = sI(vN,vM[i],vMc[i],sign);
	}

	TComplex norm,normA,suma;
	norm = (1.-v*v)/(1.-vN*vN);
	normA = normA.Power(norm,6);
	suma = normA*nor[3]*mul[0]*mul[1]*mul[2];
	return suma;
}

TComplex FFactor::GEp (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a+b+t/(4.*massP*massP)*(c+d);
}

TComplex FFactor::GMp (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a+b+c+d;
}

TComplex FFactor::GEn (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a-b+t/(4.*massP*massP)*(c-d);
}

TComplex FFactor::GMn (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a-b+c-d;
}

double FFactor::AbsGEp (TComplex t)
{
	return FFactor::GEp(t).Rho();
}

double FFactor::AbsGMp (TComplex t)
{
	return FFactor::GMp(t).Rho();
}

double FFactor::AbsGEn (TComplex t)
{
	return FFactor::GEn(t).Rho();
}

double FFactor::AbsGMn (TComplex t)
{
	return FFactor::GMn(t).Rho();
}

double FFactor::SigmaTotalP ( const TComplex &t )
{
	TComplex a = a.Abs(FFactor::GMp(t));
	TComplex b = b.Abs(FFactor::GEp(t));
	TComplex k = 4./3.*TMath::Pi()*alpha*alpha/t;
	TComplex r = r.Sqrt(1.-4*massP*massP/t);
	TComplex v =  k*r*(a*a+2.*massP*massP/t*b*b);
	return v;
}

double FFactor::SigmaTotalN ( const TComplex &t )
{
	TComplex a = a.Abs(FFactor::GMn(t));
	TComplex b = b.Abs(FFactor::GEn(t));
	TComplex k = 4./3.*TMath::Pi()*alpha*alpha/t;
	TComplex r = r.Sqrt(1.-4*massN*massN/t);
	TComplex v =  k*r*(a*a+2.*massN*massN/t*b*b);
	return v;
}

/// Options: proton, neutron, all 

void FFactor::PerformCheck ( const char* nameString, double sD )
{
	std::cout << "\n>> Standard checks for `" << nameString << "' form factors:" << std::endl;
	double chkA, chkB, chkD, chkF;
	double normaG, normaT;
	
	// Threshold
	handSome = false;

	if ( (nameString == "proton") || (nameString == "all") ) {
		chkA = FFactor::GEp(trashP).Re();
		chkB = FFactor::GMp(trashP).Re();
		chkF = shortCut::howDiff(chkA,chkB,sD);
		if ( chkF == 1 ) {
			chkD = fabs(chkA-chkB);
			std::cout << "> Warning: GEp-GMp at threshold = " << chkD << " > " << sD << std::endl;
			handSome = true;
		}
	}

	if ( (nameString == "neutron") || (nameString == "all") ) {
		chkA = FFactor::GEn(trashN).Re();
		chkB = FFactor::GMn(trashN).Re();
		chkF = shortCut::howDiff(chkA,chkB,sD);
		if ( chkF == 1 ) {
			chkD = fabs(chkA-chkB);
			std::cout << "> Warning: GEn-GMn at threshold = " << chkD << " > " << sD << std::endl;
			handSome = true;
		}
	}

	if (handSome == false) {
		std::cout << "> [OK] ... value at threshold" << std::endl;
	}
	
	// Normalization
	handSome = false;

	if ( (nameString == "proton") || (nameString == "all") ) {
		normaG = FFactor::GEp(k0).Re();
		normaT = 1.;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> Warning: GEp(0) = " << normaG << std::endl;
			handSome = true;
		}

		normaG = FFactor::GMp(k0).Re();
		normaT = 1.+ammP;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> Warning: GMp(0) = " << normaG << std::endl;
			handSome = true;
		}
	}

	if ( (nameString == "neutron") || (nameString == "all") ) {
		normaG = FFactor::GEn(k0).Re();
		normaT = 0.;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> Warning: GEn(0) = " << normaG << std::endl;
			handSome = true;
		}

		normaG = FFactor::GMn(k0).Re();
		normaT = ammN;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> Warning: GMn(0) = " << normaG << std::endl;
			handSome = true;
		}
	}

	if (handSome == false) {
		std::cout << "> [OK] ... normalization" << std::endl;
	}	
}

	}  // namespace Minuit2

}  // namespace ROOT
