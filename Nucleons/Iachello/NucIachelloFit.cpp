/**
 * $Date$
 * $Revision$
 * $Author$
 * 
 * @file
 * @brief       Fit of Iachello form factors.
 */

using namespace ROOT::Minuit2;

/// Perform fit of nucleon form factors with Minuit2.

void performFit ( char* p, char* f1, char* f2, char* o) {

	//f = "dataNucleonsApprox.dat";
	std::cout << "\n> Fit of data:" << std::endl;
	std::cout << "> Form factor data:    `" << f1 << "'" << std::endl;
	std::cout << "> Cross section data:  `" << f2 << "'" << std::endl;
	std::cout << "> Start parameters:    `" << p << "'" << std::endl;
	std::cout << "> Output parameters:   `" << o << "'" << std::endl;
	std::cout << "> Covariance matrix:   `" << matrixFile << "'" << std::endl;

	FFactor pFF(5);
	pFF.LoadParameters(p);
	int nPar = pFF.numberOfParameters;
	std::cout << "> Number of FF parameters:       " << nPar << std::endl;
	//pFF.PrintParameters();

	ExperimentalData Z;
	Z.ReadData(f1);
	Z.ReadDataCrossSection(f2);
	Z.DataInfo();
	//Z.RemoveDataType(133);
	//Z.RemoveDataType(128);
	
	int nData = Z.size();
	std::cout << "> Number of fitted points:       " << nData << std::endl;
	//Z.CheckData(); //!!!

	/// Create FCN function
	NucleonFcn fFCN(nPar, Z.Type(), Z.X(), Z.Val(), Z.ErrUp(), Z.ErrDown(), Z.Theta(), Z.Energy());

	/// Minuit
	double step = 0.01;
	MnUserParameters upar;
	for (int i = 0; i < nPar; ++i) {
		upar.Add(pFF.v[i].name, pFF.v[i].val, step, pFF.v[i].down, pFF.v[i].up);
	}

	MnStrategy(2);
	/// Create minimizer
	MnMigrad migrad(fFCN, upar);
	MnMinimize minimize(fFCN, upar);
	MnSimplex simplex(fFCN, upar);
			
	//migrad.Fix(0.);
	//migrad.Fix(1.);
	//migrad.Fix(2.);
	//migrad.Fix(3.);
	//migrad.Fix(4.);

	/// Fitting
	TStopwatch timer;
	timer.Start();

	std::cout << "\n> Minimalization:" << std::endl;

	// operator(maxfcn,tolerance):
	// maxfcn: maximum function calls
	// migrad EDM: 0.001*tolerance*up
	// simplex:    tolerance*up
	
	std::cout << "> 1. minimalization..." << std::endl;	
	//FunctionMinimum min = migrad(100000,100.);
	FunctionMinimum min = minimize(1000000,100.);

	//minimize.Fix(0.);
	//minimize.Fix(1.);
	//minimize.Fix(2.);
	//minimize.RemoveLimits(1.);
	std::cout << "> 2. minimalization..." << std::endl;
	FunctionMinimum min2 = migrad(10000,1.);

	std::cout << "> in progress..." << std::endl;
	//FunctionMinimum min9 = migrad(1000000,50.);
	FunctionMinimum min9 = minimize(1000000,50.);
	//FunctionMinimum min9 = simplex(1000000,1.);

	timer.Stop();

	/// Standard output
	std::cout << "> Minimum: " << min9 << std::endl;
	std::cout << "> FCN value:     " << min9.Fval() << std::endl;
	std::cout << "> FCN value/ndf: " << min9.Fval()/(nData-nPar) << std::endl;
	std::cout << "> Points       : " << nData << std::endl;

	std::cout << min9.UserState() << std::endl;
	std::cout << min9.UserCovariance() << std::endl;
	
	std::cout << "> Real time of minimalization: " << timer.RealTime() << " sec." << std::endl;

	/// Output parameters
	std::ofstream os;
	os.open(o);
	os.precision(8);
	for (int i = 0; i < nPar; ++i) {
		os << i+1 << ", " << min9.UserParameters().Name(i) << ", ";
		os.width(12);
		os << min9.UserParameters().Value(i) << ", ";
		os.width(12);
		os << min9.UserParameters().Error(i) << ", ";
		os << pFF.v[i].down << ", " << pFF.v[i].up << std::endl;
	}
	os.close();

	/// Covariance matrix
	MnUserCovariance cov = min9.UserCovariance();

	os.open(matrixFile);
	os.precision(10);
	for (int i = 0; i < nPar; ++i) {
		for (int j = 0; j < nPar; ++j) {
			os.width(17); os << cov(i,j);
		}
		os << "\n";
	}
	os.close();

}
