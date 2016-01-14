/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Calculate proton radius value, v. with cross section data.
 */

using namespace ROOT::Minuit2;

/**
 * Calculate final value and uncertainty of proton radius.
 * 
 * Input parameters: covariance matrix, form factor parameters.
 */

void performRadius ( char* p, char* c ) {

	std::cout << "\n> Radius calculation:   '" << std::endl;
	std::cout << "> Radius parameters:   '" << p << "'" << std::endl;
	std::cout << "> Covarince matrix:    '" << c << "'" << std::endl;	
	std::cout << "> Radius data output:  '" << radiusFile << "'\n" << std::endl;	

	FFactor Rp(12);
	Rp.LoadParameters(p);
	Rp.CheckParameters();
	Rp.LoadCovMatrix(c);
	Rp.CheckCovMatrix();
	ExperimentalData Z;
	Z.ReadCovMatrix(c,12);
	int nPar = Rp.numberOfParameters;
	
	std::cout << "\n> Covariance matrix:" << std::endl;
/* 	for (int i = 0; i < nPar; ++i) {
 * 		for (int j = 0; j < nPar; ++j) {
 * 			std::cout.precision(6);
 * 			std::cout.width(13);
 * 			std::cout << Z.cov[i][j];
 * 		}
 * 		std::cout << std::endl;
 * 	}
 */
	
	std::cout << "> Diagonal elements and parameter errors:" << std::endl;
	std::cout << " || No. | Difference | Diagonal element of cov. matrix | Parameter's error ||" << std::endl;
	for (int i = 0; i < nPar; ++i) {
		std::cout << i+1 << ". " << Z.Cov()(i,i) - Rp.E(i)*Rp.E(i) << "   " << Z.Cov()(i,i) << "   " << Rp.E(i) << std::endl;
	}
	//Rp.PrintParameters();

	std::cout << "\n || Step | Radius value     | Radius uncertainty     ||" << std::endl;
	std::cout.precision(8);
	double radVal[4], radErr[4];
	double r, radius, radius_err;
	for (int t = 0; t < 4; ++t) {
		std::cout << "Form factor type = " << t << std::endl;
		r = 0.1;
		for (int i = 0; i < 6; ++i) {
			double s = 1.;
			for (int j = 0; j < i+1; ++j) {
				s = s*r;
			}
			radius = Rp.Radius(t, s);
			radius_err = Rp.RadiusUncer(t, s);
			std::cout.width(10); std::cout <<  s;
			std::cout.width(17); std::cout << radius;
			std::cout.width(17); std::cout << radius_err << std::endl;
		}
		radVal[t] = radius;
		radErr[t] = radius_err;
	}
	//std::cout << "> Do it:" << std::endl;	
	//std::vector< std::vector<double> > em = Z.cov;
	//for (int i = 0; i < nPar; ++i) {
	//	std::cout << i+1 << ". " << em[i][i] << "   " << Z.cov[i][i] << "   " << std::endl;
	//}
	
	std::cout << "\n> Radii + MC values: " << std::endl;

	for (int t = 0; t < 4; ++t) {
		Rp.RadiusUncerMC(t, 0.0001, 10000);
		std::cout << "Cov. matrix value : ";
		std::cout.width(17); std::cout << radVal[t] << "  +/-  ";
		std::cout.width(12); std::cout << radErr[t] << std::endl;
	}

/*  	DataGenerator mydata(123), mydata2(123);
 * 	mydata.Gauss(1.,1.);
 * 	mydata2.Flat(1.,1.);
 * 	std::vector<double> value = mydata.Data();
 * 	std::vector<double> value2 = mydata2.Data();
 * 	std::cout << "\n" << std::endl;
 * 	std::cout << mydata.Size() << std::endl;
 *  	for (int i = 0; i < mydata.Size(); ++i) {
 *  		std::cout << i+1 << " " << value[i] << " " << value2[i] << std::endl;
 *  	}
 */

	
}
