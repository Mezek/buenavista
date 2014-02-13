/**
 * $Date$
 * $Revision$
 * $Author$
 * 
 * @file
 * @brief	Chi2 calculation for Iachello model.
 */

using namespace ROOT::Minuit2;

/// Perform Chi2 calculation with provided data.

void performChi ( char* p, char* f1, char* f2 ) {

	std::cout << "\n> Chi2 calculation:" << std::endl;
	std::cout << "> Chi2 parameters:         `" << p << "'" << std::endl;
	std::cout << "> Chi2 form factor data:   `" << f1 << "'" << std::endl;
	std::cout << "> Chi2 cross section data: `" << f2 << "'" << std::endl;

	ExperimentalData C;
	C.ReadData(f1);
	C.ReadDataCrossSection(f2);
	C.DataInfo();

	//C.RemoveDataType(133);
	//C.RemoveDataType(128);
	int nDat = C.size();
	
	std::vector<int> type = C.Type();
	std::vector<double> x = C.X();
	std::vector<double> val = C.Val();
	std::vector<double> errUp = C.ErrUp();
	std::vector<double> errDown = C.ErrDown();
	std::vector<double> theta = C.Theta();
	std::vector<double> energy = C.Energy();	
	
	ofstream os (xiFile);
	FFactor bornTBW(5);
	bornTBW.LoadParameters(p);
	bornTBW.CheckParameters();
	bornTBW.PrintParameters();
	int nPar = bornTBW.numberOfParameters;

	/// Alternative to calculate Xi squared

	std::vector<double> chipar(nPar);
	for (int i = 0; i < nPar; ++i) { chipar[i] = bornTBW.A(i); }
	NucleonFcn *myFunc= new NucleonFcn(nPar, type, x, val, errUp, errDown, theta, energy);
	const FCNBase& myF(*myFunc);
	double Chi2Alt = myF(chipar);
	
	double limChi = 12.;
	double chi = 0.;
	double chi2 = 0.;
	double Chi2 = 0.;
	double delta = 0.;
	double sigmaAverage = 0.;
	double Asym = 0.;
	TComplex nFv;
	double tau,eps,Es,Mott,DipFF,InvDipFF;
	
	os << "Partial Chi > " << limChi << std::endl;
	os << "n \t type[n] \t x[n] \t val[n] \t nFv \t chi2 " << std::endl;

	for(int n = 0; n < nDat; ++n) {
		nFv = bornTBW.TypeDefVal(type[n],x[n],theta[n],energy[n]);
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
			os.width(5);
			os << n+1;
			os.width(4);
			os << type[n];
			os.width(12);
			os << x[n];
			os.width(12);
			os << val[n];
			os.width(12);
			os << nFv.Re();
			os.width(12);
			os << chi2 << std::endl;
		}
	}
	os << "Chi2: " << Chi2 << std::endl;
	os.close();
	std::cout << "\n> Chi2 results:" << std::endl;
	std::cout << "Chi2    : " << Chi2 << std::endl;
	std::cout << "Points  : " << nDat << std::endl;
	std::cout << "Chi2/ndf: " << Chi2/(nDat-nPar) << std::endl;
	std::cout << "Chi2Alt : " << Chi2Alt << std::endl;
	//std::cout << cos(35.5/180.*TMath::Pi()) << " " << cos(35.5) << " " << cos(theta[1]) << std::endl;
}
