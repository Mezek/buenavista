/**
 * $Date$
 * $Revision$
 * $Author$
 * 
 * @file
 * @brief	Tabulate results for form factors in choosen interval of t.
 */

using namespace ROOT::Minuit2;

/// Perform tabulation with provided data.

void performTable ( char* p, char* f ) {

	std::cout << "\n> Tabulation:" << std::endl;
	std::cout << "> Table parameters:        `" << p << "'" << std::endl;
	std::cout << "> Table output file:       `" << f << "'" << std::endl;
	
	ofstream os (tableFile);
	FFactor bornTBW(5);
	bornTBW.LoadParameters(p);

	/// Read interval

	double tMin, tMax, tN;
	std::cout << "\n>> Set the interval of the tabulation." << std::endl;

	do {
		std::cout << ">> Variable t must be negative, t < 0" << std::endl;
		std::cout << ">> Starting point in GeV:  " << std::endl;
		std::cin >> tMin;
		std::cout << ">> End point in GeV:       " << std::endl;
		std::cin >> tMax;
	} while ((tMin > 0.) || (tMax > 0.));
	
	std::cout << ">> Number of steps:        " << std::endl;
	std::cin >> tN;

	double tStep;
	tStep = (tMax-tMin)/tN;
	double tI = 0.;

	for (int i = 0; i <= tN; i++) {
		tI = tMin + i*tStep;
		os << tI << "\t" << bornTBW.GEP(tI) << "\t" << bornTBW.GMP(tI) << "\t"
			<< bornTBW.GEN(tI) << "\t" << bornTBW.GMN(tI) << "\t" << std::endl;
	}
		
	os.close();
}
