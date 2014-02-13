/**
 * $Date$
 * $Revision$
 * $Author$
 * 
 * @file
 * @brief       Fit of pion form factors.
 */

using namespace ROOT::Minuit2;

/// Perform fit of pion form factors with Minuit2.

void performFit ( char* p, char* f, char* o) {

	std::cout << "\n> Fit of data:" << std::endl;
	std::cout << "> Form factor data:   `" << f << "'" << std::endl;
	std::cout << "> Start parameters:   `" << p << "'" << std::endl;
	std::cout << "> Output parameters:  `" << o << "'" << std::endl;

	FFactor pFF(11);
	pFF.LoadParameters(p);
	int nPar = pFF.numberOfParameters;
	std::cout << "Number of FF parameters:       " << nPar << std::endl;
	//pFF.PrintParameters();

	ExperimentalData Z;
	Z.ReadData(f);
	int nData = Z.size();
	std::cout << "Number of fitted points:       " << nData << std::endl;

	/// Create FCN function
	PionFcn fFCN(nPar, Z.Type(), Z.X(), Z.Val(), Z.ErrUp(), Z.ErrDown());

	/// Minuit
	double step = 0.01;
	MnUserParameters upar;
	for (int i = 0; i < nPar; ++i) {
		upar.Add(pFF.AName(i), pFF.AVal(i), step, pFF.ADown(i), pFF.AUp(i));
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
	/*migrad.Fix(4.);
	migrad.Fix(5.);
	migrad.Fix(6.);
	migrad.Fix(7.);
	migrad.Fix(8.);
	migrad.Fix(9.);
	migrad.Fix(10.);
	migrad.Fix(11.);
	migrad.Fix(12.);
	migrad.Fix(13.);
	migrad.Fix(14.);
	migrad.Fix(15.);*/

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

	//minimize.Fix(0.);
	//minimize.Fix(1.);
	//std::cout << "> ... 2. minimalization..." << std::endl;
	//FunctionMinimum min2 = migrad(10000,1.);

	std::cout << "> ... in progress..." << std::endl;
	//FunctionMinimum min9 = migrad(1000000,50.);
	FunctionMinimum min9 = minimize(1000000,50.);
	//FunctionMinimum min9 = simplex(1000000,1.);

	timer.Stop();

	/// Standard output
	std::cout << "> Minimum: " << min9 << std::endl;
	std::cout << "> FCN value: " << min9.Fval() << std::endl;
	std::cout << "> Real time of minimalization: " << timer.RealTime() << " sec." << std::endl;

	//pFF.PrintParameters();
	
	//~ // Scan and Plot parameters
 	//~ {
    //~ MnScan scan(fFCN, upar, 1);
    //~ std::cout << "Scan parameters: " << scan.Parameters() << std::endl;
    //~ MnPlot plot;
    //~ for (unsigned int i = 0; i < upar.VariableParameters(); ++i) {
		//~ std::vector<std::pair<double, double> > xy = scan.Scan(i);
		//~ //std::vector<std::pair<double, double> > xy = scan.scan(0);
		//~ std::cout << i+1 << ". parameter - `" << upar.Name(i) << "'"<< std::endl;
		//~ plot(xy);
	//~ }
    //~ std::cout << scan.Parameters() << std::endl;
	//~ }

	// Output parameters
	ofstream myOutputParam (o);
	for (int i = 0; i < nPar; ++i) {
		myOutputParam << i+1 << ", " << min9.UserParameters().Name(i) << ", ";
		myOutputParam << std::setprecision(10) << min9.UserParameters().Value(i) << ", ";
		myOutputParam << std::setprecision(10) << min9.UserParameters().Error(i) << ", ";
		myOutputParam << pFF.ADown(i) << ", " << pFF.AUp(i) << std::endl;
	}
	myOutputParam.close();

}
