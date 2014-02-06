/**
 * $Date: 2013-06-19 12:42:41 +0200 (Wed, 19 Jun 2013) $
 * $Revision: 377 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/NucleonsChiH.cpp $
 * $Id: NucleonsChiH.cpp 377 2013-06-19 10:42:41Z bartos $
 *
 * @file
 * @brief	Chi2 calculation for nucleon form factors, v. hyperon.
 */

using namespace ROOT::Minuit2;

/// Perform Chi2 calculation with provided data.

void performChi ( char* p, char* f ) {

	std::cout << "\n> Chi2 calculation:" << std::endl;
	std::cout << "> Chi2 parameters:    `" << p << "'" << std::endl;
	std::cout << "> Chi2 data:          `" << f << "'" << std::endl;

	ExperimentalData C;
	C.ReadData(f);
	C.DataInfo();

	int nDat = C.size();
	
	std::vector<int> type = C.Type();
	std::vector<double> x = C.X();
	std::vector<double> val = C.Val();
	std::vector<double> errUp = C.ErrUp();
	std::vector<double> errDown = C.ErrDown();
	std::vector<double> theta = C.Theta();
	std::vector<double> energy = C.Energy();	
	
	ofstream os (xiFile);
	FFactor bornTBW;
	bornTBW.LoadParameters(p);
	bornTBW.CheckParameters();
	bornTBW.PrintParameters();
	int nPar = bornTBW.numberOfParameters;

	double limChi = 12.;
	double chi = 0.;
	double chi2 = 0.;
	double Chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;
	double Asym = 0.;
	TComplex nFv,nE,nM;
	os << "Partial Chi > " << limChi << std::endl;
	for (int n = 0; n < nDat; ++n) {
		// "protonElectric"
		if ((type[n] == 1) || (type[n] == 2)) {
			nFv = bornTBW.AbsGE1(x[n]);
			//nFv = bornTBW.ScalarOne(x[n]);
			//nFv = bornTBW.ScalarTwo(x[n]);
			//nFv = bornTBW.VectorOne(x[n]);
			//nFv = bornTBW.VectorTwo(x[n]);
		}
		// "protonMagnetic"
		if ((type[n] == 3) || (type[n] == 4)) {
			nFv = bornTBW.AbsGM1(x[n]);
		}
		// "neutronElectric"
		if ((type[n] == 5) || (type[n] == 6)) {
			nFv = bornTBW.AbsGE2(x[n]);
		}
		// "neutronMagnetic"
		if (type[n] == 7) {
			nFv = bornTBW.GM2(x[n]);
		}
		if (type[n] == 8) {
			nFv = bornTBW.AbsGM2(x[n]);
		}
		// "protonRatios"
		if (type[n] == 9) {
			nE = bornTBW.GE1(x[n]);
			nM = bornTBW.GM1(x[n]);
			nFv = nFv.Abs(muP*nE/nM);
		}
		// "neutronRatios"
		if (type[n] == 10) {
			nE = bornTBW.GE2(x[n]);
			nM = bornTBW.GM2(x[n]);
			nFv = nFv.Abs(muN*nE/nM);
		}

		delta = nFv.Re()-val[n];
		sigmaAverage = (errUp[n] + errDown[n])/2.;
		Asym = (errUp[n] - errDown[n])/(errUp[n] + errDown[n]);
		chi = delta/sigmaAverage;
		if (Asym == 0.) { chi2 = chi*chi; } else {
			//chi2 = chi*chi-2.*Asym*chi*chi*chi+5.*Asym*Asym*chi*chi*chi*chi;
			chi2 = chi*chi;
		}
		Chi2 += chi2;
		if ( chi2 > limChi ) {
			os << "\n" << n << "\t" << type[n] << "\t" << x[n] << "\t" << val[n] << "\t" << chi2 << std::endl;
		}
	}
	os << "Chi2: " << Chi2 << std::endl;
	os.close();
	std::cout << "\n> Chi2 results:" << std::endl;
	std::cout << "Chi2    : " << Chi2 << std::endl;
	std::cout << "Points  : " << nDat << std::endl;
	std::cout << "Chi2/ndf: " << Chi2/(nDat-nPar) << std::endl;
}
