/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/Nucleon3G.cpp $
 * $Id$
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
		
FFactor::FFactor ( std::size_t size ): a(size), v(size)
{ 
	t(1.,0.00001);
	//std::cout << ">> FFactor is empty! Value of t: " << t << std::endl;

	// izo :: name, nor, mesons

	// Names
	FF[0].name = "F1s";
	FF[1].name = "F1v";
	FF[2].name = "F2s";
	FF[3].name = "F2v";
		
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
	std::cout << "F2v\t" << "\t" << "\t"<< "\t" << "\t" << "\t" << std::endl;

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
	//std::cout << "F1s: v = " << v << " vN = " << vN << std::cout << std::endl;
	
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

TComplex FFactor::GEP (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a+b+t/(4.*massP*massP)*(c+d);
}

TComplex FFactor::GMP (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a+b+c+d;
}

TComplex FFactor::GEN (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a-b+t/(4.*massN*massN)*(c-d);
}

TComplex FFactor::GMN (TComplex t)
{
	TComplex a = FFactor::ScalarOne(t);
	TComplex b = FFactor::VectorOne(t);
	TComplex c = FFactor::ScalarTwo(t);
	TComplex d = FFactor::VectorTwo(t);
	return a-b+c-d;
}

double FFactor::AbsGEP (TComplex t)
{
	return FFactor::GEP(t).Rho();
}

double FFactor::AbsGMP (TComplex t)
{
	return FFactor::GMP(t).Rho();
}

double FFactor::AbsGEN (TComplex t)
{
	return FFactor::GEN(t).Rho();
}

double FFactor::AbsGMN (TComplex t)
{
	return FFactor::GMN(t).Rho();
}

double FFactor::SigmaTotalP ( const TComplex &t )
{
	TComplex a = a.Abs(FFactor::GMP(t));
	TComplex b = b.Abs(FFactor::GEP(t));
	TComplex k = 4./3.*TMath::Pi()*alpha*alpha/t;
	TComplex r = r.Sqrt(1.-4*massP*massP/t);
	TComplex v = k*r*(a*a+2.*massP*massP/t*b*b);
	return v;
}

double FFactor::SigmaTotalN ( const TComplex &t )
{
	TComplex a = a.Abs(FFactor::GMN(t));
	TComplex b = b.Abs(FFactor::GEN(t));
	TComplex k = 4./3.*TMath::Pi()*alpha*alpha/t;
	TComplex r = r.Sqrt(1.-4*massN*massN/t);
	TComplex v = k*r*(a*a+2.*massN*massN/t*b*b);
	return v;
}

TComplex FFactor::TypeDefVal ( const int &type, const TComplex &t, const double &theta, const double &energy )
{
	TComplex val,nE,nM;
	double tau,eps,th,Es,Mott,DipFF,InvDipFF;

	/** "protonElectric" */
	if ((type == 1) || (type == 2)) { val = FFactor::AbsGEP(t);	}
	/** "protonMagnetic" */
	if ((type == 3) || (type == 4)) { val = FFactor::AbsGMP(t);	}
	/** "neutronElectric" */
	if ((type == 5) || (type == 6)) { val = FFactor::AbsGEN(t);	}
	/** "neutronMagnetic" */
	if (type == 7) { val = FFactor::GMN(t); }
	if (type == 8) { val = FFactor::AbsGMN(t); }
	/** "protonRatios" */
	if (type == 9) {
		nE = FFactor::AbsGEP(t);
		nM = FFactor::AbsGMP(t);
		val = val.Abs((1.+ammP)*nE/nM);
	}
	/** "neutronRatios" */
	if (type == 10) {
		nE = FFactor::AbsGEN(t);
		nM = FFactor::AbsGMN(t);
		val = val.Abs(ammN*nE/nM);
	}
	/** "Mainz Ratios" */
	if (type >= 100) {		
		th = theta/180.*TMath::Pi();
		nE = FFactor::AbsGEP(t);
		nM = FFactor::AbsGMP(t);
		tau = t/(4*massP*massP);
		eps = 1./(1.+2.*(1.+tau)*tan(th/2.)*tan(th/2.));
		Es = energy/(1.+energy/massP*(1.-cos(th)));
		InvDipFF = (1.-t/0.71)*(1.-t/0.71);
		DipFF = 1./InvDipFF;
		val = (eps*nE*nE+tau*nM*nM)/(eps+tau*(1.+ammP)*(1.+ammP))*InvDipFF*InvDipFF;
	}
	return val;
}

/// First derivation in point t, with arbitrary step

double FFactor::Derive ( const int fftype, const TComplex &t, const double step )
{
	TComplex zp,zm;
	switch (fftype) {
		case 1:
			zp = FFactor::GMP(t+step);
			zm = FFactor::GMP(t-step);
			break;
		case 2:
			zp = FFactor::GEN(t+step);
			zm = FFactor::GEN(t-step);
			break;
		case 3:
			zp = FFactor::GMN(t+step);
			zm = FFactor::GMN(t-step);
			break;
		default:
			zp = FFactor::GEP(t+step);
			zm = FFactor::GEP(t-step);
	}
	TComplex r = (zp-zm)/(2.*step);
	return r;
}

double FFactor::DeriveOld ( const TComplex &t, const double step )
{
	TComplex zp,zm;
	zp = FFactor::GEP(t+step);
	zm = FFactor::GEP(t-step);
	TComplex r = (zp-zm)/(2.*step);
	return r;
}

/// Second derivation according to t and a_i parameter, in point t, with arbitrary step

double FFactor::DeriveXA ( const int fftype, const TComplex &t, int i, const double step )
{
	double b;
	TComplex zp,zm;
	b = a[i].val;
	FFactor::SetParameter(i,b+step);
	zp = FFactor::GEP(t+step) - FFactor::GEP(t-step);
	FFactor::SetParameter(i,b-step);
	zm = FFactor::GEP(t+step) - FFactor::GEP(t-step);
	FFactor::SetParameter(i,b);
	TComplex r = (zp-zm)/(4.*step*step);
	return r;
}

double FFactor::DeriveXAOld ( const TComplex &t, int i, const double step )
{
	double b;
	TComplex zp,zm;
	b = a[i].val;
	FFactor::SetParameter(i,b+step);
	zp = FFactor::GEP(t+step) - FFactor::GEP(t-step);
	FFactor::SetParameter(i,b-step);
	zm = FFactor::GEP(t+step) - FFactor::GEP(t-step);
	FFactor::SetParameter(i,b);
	TComplex r = (zp-zm)/(4.*step*step);
	return r;
}

double FFactor::RadiusEP ( const double step )
{
	if (step > t0v) {
		std::cout << "\n> RadiusEP: Warning! Step = " << step << " must be lower than t0v = " << t0v << std::endl;
	}
	TComplex d = FFactor::DeriveOld(0.,step);
	TComplex v2 = 6.*d*0.1*hTransC2;
	double v = sqrt(v2.Re());
	return v;
}

double FFactor::RadiusEPUncer ( const double step )
{
	// Check step value
	if (step > t0v) {
		std::cout << "\n> RadiusEPUncer: Warning! Step = " << step << " must be lower than t0v = " << t0v << std::endl;
	}

	double d[modelPar];
	TVectorD vecE(modelPar);

	for (int i = 0; i < modelPar; ++i) {
		d[i] = fabs(FFactor::DeriveXAOld(0.,i,step));
		vecE(i) = d[i];
		//std::cout << d[i] << " ";
	}
	//std::cout << std::endl;

	double sig = 0.;
	for (int i = 0; i < modelPar; ++i) {
		for (int j = 0; j < modelPar; ++j) {
			sig = sig + d[i]*d[j]*cov(i,j);
		}
	}

	// Alternative calculation only with vectors
	TVectorD vecT = cov*vecE;
	double sig2 = vecE*vecT;

	/*double sig = 0.;
	for (int i = 0; i < modelPar; ++i) {
		sig = sig + d[i]*d[i]*cov(i,i);
	}*/

	//double v = sqrt(sig)*3.*0.1*hTransC2/FFactor::RadiusEP(step);
	double v = sqrt(sig)*2./FFactor::RadiusEP(step);
	return v;
}

void FFactor::RadiusEPUncerMC ( const double step, const double nIter )
{
	// Check step value
	if (step > t0v) {
		std::cout << "\n> RadiusEPUncerMC: Warning! Step = " << step << " must be lower than t0v = " << t0v << std::endl;
	}

	// Get Cholesky decomposition of correlation matrix
	TDecompChol matA(cov);
	matA.Decompose();

	// Get lower matrix
	TMatrixD matB(modelPar,modelPar);
	matB=matA.GetU();
	matB.Transpose(matB);
	//matB.Print();

	TVectorD parMean(modelPar), parSigma(modelPar);

	for (int i = 0; i < modelPar; ++i) {
		parMean(i) = FFactor::A(i);
		// Sigma of parameters, transform too?!
		parSigma(i) = FFactor::E(i)*FFactor::E(i); // parSigma(i) = Z.cov[i][i];
	}

	TRandom3 gen_val;
	//gen_val.SetSeed(0);
	//std::cout << gen_val.GetSeed() << std::endl;
	TVectorD vecX(modelPar);
	double mc_mean = FFactor::RadiusEP(step);
	double mc_sigma;
	double mc_variance = 0.;
	double mc_mean_check = 0.;
	int mc_num = nIter;
	for (int k = 0; k < mc_num; ++k) {
		std::vector<double> tempVal, tempValX;
		tempVal.resize(modelPar);
		tempValX.resize(modelPar);
		for (int i = 0; i < modelPar; ++i) {
			//vecX[i] = gen_val.Gaus(parMean(i),parSigma(i));
			vecX(i) = gen_val.Gaus(0.,1.);
			tempValX[i] = vecX(i);
		}
		// Transform vector by matrix
		vecX *= matB; // implemented: vector = matrix*vector
		// (Un)correlated values
		for (int i = 0; i < modelPar; ++i) {
			tempVal[i] = vecX(i)+parMean(i);
			//std::cout << tempValX[i] << " : " << vecX(i) << " : " << parMean(i) << " : " << tempVal[i] << std::endl;
		}
		FFactor::SetParameters(tempVal);
		double x_i = FFactor::RadiusEP(step);
		if (k%1000 == 0) {
			std::cout << "\r>> " << k/1000+1 << " / " << mc_num/1000 << "k";
			//std::cout << "= \b";
			fflush(stdout);
		}
		mc_variance = (x_i-mc_mean)*(x_i-mc_mean) + mc_variance;
		mc_mean_check = mc_mean_check + x_i;
	}
	mc_sigma = sqrt(mc_variance/(mc_num-1));
	mc_mean_check = mc_mean_check/mc_num;

	std::cout << "\n>> Proton radius value : Proton radius uncertainty" << std::endl;
	std::cout.width(17); std::cout << mc_mean;
	std::cout.width(17); std::cout << sqrt(mc_sigma) << " " << mc_sigma/sqrt(mc_num) << std::endl;
	std::cout << ">> Proton radius value check" << std::endl;
	std::cout.width(17); std::cout << mc_mean_check;
}


/// Options: proton, neutron, all

void FFactor::CheckFormFactor ( const char* nameString, double sD )
{
	std::cout << "\n>> Standard checks for `" << nameString << "' form factors:" << std::endl;
	double chkA, chkB, chkD, chkF;
	double normaG, normaT;
	
	// Threshold
	handSome = false;

	if ( (nameString == "proton") || (nameString == "all") ) {
		chkA = FFactor::GEP(trashP).Re();
		chkB = FFactor::GMP(trashP).Re();
		chkF = shortCut::howDiff(chkA,chkB,sD);
		if ( chkF == 1 ) {
			chkD = fabs(chkA-chkB);
			std::cout << "> CheckFormFactor: Warning!" << std::endl;
			std::cout << "> GEp-GMp at threshold = " << chkD << " > " << sD << std::endl;
			handSome = true;
		}
	}

	if ( (nameString == "neutron") || (nameString == "all") ) {
		chkA = FFactor::GEN(trashN).Re();
		chkB = FFactor::GMN(trashN).Re();
		chkF = shortCut::howDiff(chkA,chkB,sD);
		if ( chkF == 1 ) {
			chkD = fabs(chkA-chkB);
			std::cout << "> CheckFormFactor: Warning!" << std::endl;
			std::cout << "> GEn-GMn at threshold = " << chkD << " > " << sD << std::endl;
			handSome = true;
		}
	}

	if (handSome == false) {
		std::cout << "> [OK] ... value at threshold" << std::endl;
	}
	
	// Normalization
	handSome = false;

	if ( (nameString == "proton") || (nameString == "all") ) {
		normaG = FFactor::GEP(k0).Re();
		normaT = 1.;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> CheckFormFactor: Warning! GEp(0) = " << normaG << std::endl;
			handSome = true;
		}

		normaG = FFactor::GMP(k0).Re();
		normaT = 1.+ammP;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> CheckFormFactor: Warning! GMp(0) = " << normaG << std::endl;
			handSome = true;
		}
	}

	if ( (nameString == "neutron") || (nameString == "all") ) {
		normaG = FFactor::GEN(k0).Re();
		normaT = 0.;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> CheckFormFactor: Warning! GEn(0) = " << normaG << std::endl;
			handSome = true;
		}

		normaG = FFactor::GMN(k0).Re();
		normaT = ammN;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> CheckFormFactor: Warning! GMn(0) = " << normaG << std::endl;
			handSome = true;
		}
	}

	if (handSome == false) {
		std::cout << "> [OK] ... normalization" << std::endl;
	}	
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
