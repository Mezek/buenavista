/**
 * $Date$
 * $Revision$
 * $Author$
 * 
 * @file
 * @brief	Chi2 calculation for pion form factors.
 */

using namespace ROOT::Minuit2;

/// Perform Chi2 calculation with provided data.

void performChi ( char* p, char* f ) {

	std::cout << "\n> Chi2 calculation:" << std::endl;
	std::cout << "> Chi2 parameters:    `" << p << "'" << std::endl;
	std::cout << "> Chi2 data:          `" << f << "'" << std::endl;

	ExperimentalData C;
	C.ReadData(f);
	std::vector<int> type = C.Type();
	std::vector<double> x = C.X();
	std::vector<double> val = C.Val();
	std::vector<double> errUp = C.ErrUp();
	std::vector<double> errDown = C.ErrDown();
	
	ofstream myOutputXi (xiFile);
	FFactor xiFF(12);
	xiFF.LoadParameters(p);
	xiFF.CheckParameters();
	xiFF.PrintParameters();
	int nP = C.size();
	double limChi = 7.;
	double chi = 0.;
	double chi2 = 0.;
	double Chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;
	TComplex PV;
	TComplex pFv;
	
	myOutputXi << "Partial Chi > " << limChi << std::endl;
	for(int n = 0; n < nP; ++n) {
		/** "protonElectric" */
		if ((type[n] == 0) || (type[n] == 1)) {
			PV = xiFF.Value(x[n]);
			pFv = pFv.Abs(PV);

		}

		delta = pFv.Re()-val[n];
		sigmaAverage = (errUp[n] + errDown[n])/2.;
		chi = delta/sigmaAverage;
		chi2 = chi*chi;
		
		Chi2 += chi2;
		if ( chi2 > limChi ) {
			myOutputXi << "\n" << n << "\t" << type[n] << "\t" << x[n] << "\t" << val[n] << "\t" << chi2 << std::endl;
		}
	}
	myOutputXi << "Chi2: " << Chi2 << std::endl;
	myOutputXi.close();
	std::cout << "\n> Chi2  : " << Chi2 << std::endl;
	std::cout << "> Points: " << nP << std::endl;
}
