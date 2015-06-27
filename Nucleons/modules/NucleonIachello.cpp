/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/NucleonIachello.cpp $
 * $Id$
 *
 * @file
 * @brief	Model of Iachello form factor for nucleons.
 */

#include "NucleonIachello.h"

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

}  // namespace shortCut


namespace ROOT {

	namespace Minuit2 {
		
FFactor::FFactor ( std::size_t size ): a(size), v(size)
{ 
	t = 0.;
	//std::cout << ">> FFactor is empty! Value of t: " << t << std::endl;

	// izo :: name, nor, mesons, tresh

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

double FFactor::Gfun (double t)
{
	double res = 1./(1.-a[4].val*t)/(1.-a[4].val*t);
	return res;
}

double FFactor::Afun (double t)
{
	double termA = 2./TMath::Pi()*sqrt((t-4.*massPi2)/t);
	double termB = log((sqrt(4.*massPi2-t)+sqrt(-t))/2./massPi);
	return termA*termB;
}

double FFactor::Ttrans (double t)
{
	double nom = massRho*massRho+8.*widthRho*massPi/TMath::Pi();
	double den = massRho*massRho-t+(4.*massPi2-t)*widthRho*FFactor::Afun(t)/massPi;
	return nom/den;
}

double FFactor::ScalarOne (double t)
{
	double suma = 0.;
	double termA, termB, termC;
	termA = 1.-a[1].val-a[2].val;
	termB = a[1].val*massOm*massOm/(massOm*massOm-t);
	termC = a[2].val*massPhi*massPhi/(massPhi*massPhi-t);
	suma = 0.5*FFactor::Gfun(t)*(termA+termB+termC);
	return suma;
}

double FFactor::ScalarTwo (double t)
{
	double suma = 0.;
	double termA, termB;
	termA = (-0.12-a[3].val)*massOm*massOm/(massOm*massOm-t);
	termB = a[3].val*massPhi*massPhi/(massPhi*massPhi-t);
	suma = 0.5*FFactor::Gfun(t)*(termA+termB);
	return suma;
}

double FFactor::VectorOne (double t)
{
	double suma = 0.;
	suma = 0.5*FFactor::Gfun(t)*((1.-a[0].val)+a[0].val*FFactor::Ttrans(t));
	return suma;
}

double FFactor::VectorTwo (double t)
{
	double suma = 0.;
	suma = 0.5*FFactor::Gfun(t)*3.706*FFactor::Ttrans(t);
	return suma;
}

double FFactor::GEP (double t)
{
	double a = FFactor::ScalarOne(t);
	double b = FFactor::VectorOne(t);
	double c = FFactor::ScalarTwo(t);
	double d = FFactor::VectorTwo(t);
	return a+b+t/(4.*massP*massP)*(c+d);
}

double FFactor::GMP (double t)
{
	double a = FFactor::ScalarOne(t);
	double b = FFactor::VectorOne(t);
	double c = FFactor::ScalarTwo(t);
	double d = FFactor::VectorTwo(t);
	return a+b+c+d;
}

double FFactor::GEN (double t)
{
	double a = FFactor::ScalarOne(t);
	double b = FFactor::VectorOne(t);
	double c = FFactor::ScalarTwo(t);
	double d = FFactor::VectorTwo(t);
	return a-b+t/(4.*massP*massP)*(c-d);
}

double FFactor::GMN (double t)
{
	double a = FFactor::ScalarOne(t);
	double b = FFactor::VectorOne(t);
	double c = FFactor::ScalarTwo(t);
	double d = FFactor::VectorTwo(t);
	return a-b+c-d;
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
	if ((type == 1) || (type == 2)) { val = FFactor::GEP(t);	}
	/** "protonMagnetic" */
	if ((type == 3) || (type == 4)) { val = FFactor::GMP(t);	}
	/** "neutronElectric" */
	if ((type == 5) || (type == 6)) { val = FFactor::GEN(t);	}
	/** "neutronMagnetic" */
	if (type == 7) { val = FFactor::GMN(t); }
	if (type == 8) { val = FFactor::GMN(t); }
	/** "protonRatios" */
	if (type == 9) {
		nE = FFactor::GEP(t);
		nM = FFactor::GMP(t);
		val = val.Abs((1.+ammP)*nE/nM);
	}
	/** "neutronRatios" */
	if (type == 10) {
		nE = FFactor::GEN(t);
		nM = FFactor::GMN(t);
		val = val.Abs(ammN*nE/nM);
	}
	/** "Mainz Ratios" */
	if (type >= 100) {		
		th = theta/180.*TMath::Pi();
		nE = FFactor::GEP(t);
		nM = FFactor::GMP(t);
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

double FFactor::Derive ( const TComplex &t, const double step )
{
	TComplex zp,zm;
	zp = FFactor::GEP(t+step);
	zm = FFactor::GEP(t-step);
	TComplex r = (zp-zm)/(2.*step);
	return r;
}

/// Second derivation according to t and a_i parameter, in point t, with arbitrary step

double FFactor::DeriveXA ( const TComplex &t, int i, const double step )
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
	TComplex d = FFactor::Derive(0.,step);
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
		d[i] = fabs(FFactor::DeriveXA(0.,i,step));
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


	// Normalization
	handSome = false;

	if ( (nameString == "proton") || (nameString == "all") ) {
		normaG = FFactor::GEP(0.);
		normaT = 1.;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> CheckFormFactor: Warning! GEp(0) = " << normaG << std::endl;
			handSome = true;
		}

		normaG = FFactor::GMP(0.);
		normaT = 1.+ammP;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> CheckFormFactor: Warning! GMp(0) = " << normaG << std::endl;
			handSome = true;
		}
	}

	if ( (nameString == "neutron") || (nameString == "all") ) {
		normaG = FFactor::GEN(0.);
		normaT = 0.;
		chkF = shortCut::howDiff(normaG, normaT,sD);
		if ( chkF == 1) {
			std::cout << "> CheckFormFactor: Warning! GEn(0) = " << normaG << std::endl;
			handSome = true;
		}

		normaG = FFactor::GMN(0.);
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
