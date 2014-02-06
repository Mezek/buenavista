/**
 * $Date: 2013-05-22 09:44:15 +0200 (Wed, 22 May 2013) $
 * $Revision: 346 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/NucleonsFitMC.cpp $
 * $Id: NucleonsFitMC.cpp 346 2013-05-22 07:44:15Z bartos $
 * 
 * @file
 * @brief       Fit of nucleon form factors, v. MC.
 */

using namespace ROOT::Minuit2;

/// Perform fit of nucleon form factors with Minuit2.

void performFit ( char* p, char* f1, char* f2, char* o, int nSamp) {

	//f = "dataNucleonsApprox.dat";
	std::cout << "\n> Fitted data:        " << std::endl;
	std::cout << "> Form factor data:    '" << f1 << "'" << std::endl;
	std::cout << "> Cross section data:  '" << f2 << "'" << std::endl;
	std::cout << "> Start parameters:    '" << p << "'" << std::endl;
	std::cout << "> Output parameters:   '" << o << "'" << std::endl;

	FFactor pFF(12);
	pFF.LoadParameters(p);
	int nPar = pFF.numberOfParameters;

	ExperimentalData Z;
	Z.ReadData(f1);
	Z.ReadDataCrossSection(f2);
	int nData = Z.size();
	std::cout << "> Number of samples:             " << nSamp << std::endl;
	std::cout << "> Number of fitted points:       " << nData << std::endl;
	std::vector<int> type = Z.Type();
	std::vector<double> x = Z.X();
	std::vector<double> val = Z.Val();
	std::vector<double> errUp = Z.ErrUp();
	std::vector<double> errDown = Z.ErrDown();
	std::vector<double> theta = Z.Theta();
	std::vector<double> energy = Z.Energy();

	/// Fitting
	TStopwatch timer;
	timer.Start();

	double delta,sigma;
	double radius;
	std::vector<double> rand_val(nData);
	std::vector<double> new_par(nPar);

	/// Output parameters
	std::ofstream os;
	//os.open(o, std::ios::app); // To append data: ofstream is default "ios::out", we add "ios::app" 
	os.open(o);
	int k = 0;
	do {
		/// Generate data
		TRandom3 gen_val(0);
		for (int i = 0; i < nData; ++i) {
			delta = fabs(errDown[i] + errUp[i])/2.;
			sigma = delta*delta;
			rand_val[i] = gen_val.Gaus(val[i],sigma);
			//std::cout.width(12);
			//if ( (i==5) || (i==155)) { std::cout << rand_val[i] << " "; }
		}
		//std::cout << std::endl;
		//std::cout << k << " " << nSamp << std::endl;
		
		/// Create FCN function
		//NucleonFcn fFCN(nPar, Z.Type(), Z.X(), rand_val, Z.ErrUp(), Z.ErrDown(), Z.Theta(), Z.Energy());
		NucleonFcn *fFCN = new NucleonFcn(nPar, Z.Type(), Z.X(), rand_val, Z.ErrUp(), Z.ErrDown(), Z.Theta(), Z.Energy());

		/// Minuit
		double step = 0.01;
		MnUserParameters upar;
		for (int i = 0; i < nPar; ++i) {
			upar.Add(pFF.v[i].name, pFF.v[i].val, step, pFF.v[i].down, pFF.v[i].up);
		}

		MnStrategy(2);
		// Create minimizer
		//MnMigrad migrad(fFCN, upar);
		MnMinimize minimize(*fFCN, upar);
		//MnSimplex simplex(fFCN, upar);

		// operator(maxfcn,tolerance):
		// maxfcn: maximum function calls
		// migrad EDM: 0.001*tolerance*up
		// simplex:    tolerance*up
	
		std::cout << "> " << k+1 << ". minimalization..." << std::endl;	
		//FunctionMinimum min = migrad(100000,100.);
		//FunctionMinimum min = minimize(1000000,100.);
		//FunctionMinimum min9 = migrad(1000000,50.);
		FunctionMinimum min9 = minimize(1000000,50.);
		//FunctionMinimum min9 = simplex(1000000,1.);

		for (int i = 0; i < nPar; ++i) {
			new_par[i] = min9.UserParameters().Value(i);
			/*
			os << i+1 << ", " << min9.UserParameters().Name(i) << ", ";
			os.width(12);
			os << min9.UserParameters().Value(i) << ", ";
			os.width(12);
			os << min9.UserParameters().Error(i) << ", ";
			os << pFF.v[i].down << ", " << pFF.v[i].up << std::endl;
			*/
		}

		pFF.SetParameters(new_par);
		radius = pFF.RadiusEP(0.0001);
		std::cout << "> Proton radius: " << radius << std::endl;
		os.precision(10);
		os.width(15);
		os << radius << std::endl;
		delete fFCN;
		++k;

		timer.Stop();
	}
	while ( k < nSamp );
	os.close();
}
