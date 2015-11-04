/**
 * $Date$
 * $Revision$
 * $Author$
 * 
 * @file
 * @brief	Chi2 calculation for nucleon form factors.
 */

using namespace ROOT::Minuit2;

/// Perform Chi2 calculation with provided data.

void performChi ( char* p, char* f ) {

	std::cout << "\n> Chi2 calculation:" << std::endl;
	std::cout << "> Chi2 parameters:    `" << p << "'" << std::endl;
	std::cout << "> Chi2 data:          `" << f << "'" << std::endl;

	ExperimentalData C;
	C.ReadData(f);
	int nDat = C.size();

	std::vector<int> type = C.Type();
	std::vector<double> x = C.X();
	std::vector<double> val = C.Val();
	std::vector<double> errUp = C.ErrUp();
	std::vector<double> errDown = C.ErrDown();
	std::vector<double> theta = C.Theta();
	std::vector<double> energy = C.Energy();	
	
	ofstream myOutputXi (xiFile);
	FFactor bornTBW(12);
	bornTBW.LoadParameters(p);
	bornTBW.CheckParameters();
	bornTBW.PrintParameters();
	int nPar = bornTBW.numberOfParameters;

	double limChi = 50.;
	double chi = 0.;
	double chi2 = 0.;
	double Chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;
	double Asym = 0.;
	TComplex nFv,nE,nM;
	myOutputXi << "Partial Chi > " << limChi << std::endl;
	for(int n = 0; n < nDat; ++n) {
		/** "protonElectric" */
		if ((type[n] == 1) || (type[n] == 2)) {
			nFv = bornTBW.AbsGEP(x[n]);
			//nFv = bornTBW.ScalarOne(x[n]);
			//nFv = bornTBW.ScalarTwo(x[n]);
			//nFv = bornTBW.VectorOne(x[n]);
			//nFv = bornTBW.VectorTwo(x[n]);
		}
		/** "protonMagnetic" */
		if ((type[n] == 3) || (type[n] == 4)) {
			nFv = bornTBW.AbsGMP(x[n]);
		}
		/** "neutronElectric" */
		if ((type[n] == 5) || (type[n] == 6)) {
			nFv = bornTBW.AbsGEN(x[n]);
		}
		/** "neutronMagnetic" */
		if (type[n] == 7) {
			nFv = bornTBW.GMN(x[n]);
		}
		if (type[n] == 8) {
			nFv = bornTBW.AbsGMN(x[n]);
		}
		/** "protonRatios" */
		if (type[n] == 9) {
			nE = bornTBW.GEP(x[n]);
			nM = bornTBW.GMP(x[n]);
			nFv = nFv.Abs((1.+ammP)*nE/nM);
		}
		/** "neutronRatios" */
		if (type[n] == 10) {
			nE = bornTBW.GEN(x[n]);
			nM = bornTBW.GMN(x[n]);
			nFv = nFv.Abs(ammN*nE/nM);
		}

		delta = nFv.Re()-val[n];
		sigmaAverage = (errUp[n] + errDown[n])/2.;
		chi = delta/sigmaAverage;
		chi2 = chi*chi;

		/*
		delta = nFv.Re()-val[n];
		sigmaAverage = (errUp[n] + errDown[n])/2.;
		Asym = (errUp[n] - errDown[n])/(errUp[n] + errDown[n]);
		chi = delta/sigmaAverage;
		if (Asym == 0.) { chi2 = chi*chi; } else {
			chi2 = chi*chi-2.*Asym*chi*chi*chi+5.*Asym*Asym*chi*chi*chi*chi;
			//chi2 = chi*chi;
		}
		*/
		
		Chi2 += chi2;
		if ( chi2 > limChi ) {
			myOutputXi << "\n" << n << "\t" << type[n] << "\t" << x[n] << "\t" << val[n] << "\t" << chi2 << std::endl;
		}
		std::cout << std::setw(6) << std::right << n
                  << std::setw(6) << std::right << type[n]
				  << std::setw(12) << std::right << x[n]
				  << std::setw(12) << std::right << val[n]
				  << std::setw(12) << std::right << std::setprecision(4) << std::fixed << chi2
				  << std::endl;
	}
	myOutputXi << "Chi2: " << Chi2 << std::endl;
	myOutputXi.close();
	std::cout << "\n> Chi2  : " << Chi2 << std::endl;
	std::cout << "> Points  : " << nDat << std::endl;
	std::cout << "> Chi2/ndf: " << Chi2/(nDat-nPar) << std::endl;
}
