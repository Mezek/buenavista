/**
 * $Date: 2013-10-02 13:45:42 +0200 (Wed, 02 Oct 2013) $
 * $Revision: 383 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/NucleonsTableCS.cpp $
 * $Id: NucleonsTableCS.cpp 383 2013-10-02 11:45:42Z bartos $
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
	FFactor bornTBW(12);
	bornTBW.LoadParameters(p);

	/// Read interval

	double tMin, tMax, tN;
	std::cout << "\n>> Set the interval of the tabulation:" << std::endl;
	std::cout << ">> Starting point in GeV:  " << std::endl;
	std::cin >> tMin;
	std::cout << ">> End point in GeV:       " << std::endl;
	std::cin >> tMax;
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
