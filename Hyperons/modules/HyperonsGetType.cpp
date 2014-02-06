/**
 * $Date: 2013-05-21 10:51:25 +0200 (Tue, 21 May 2013) $
 * $Revision: 340 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/modules/HyperonsGetType.cpp $
 * $Id: HyperonsGetType.cpp 340 2013-05-21 08:51:25Z bartos $
 *
 * @file
 * @brief	Get type of a hyperon.
 */

const int opt = 5;

int HyperonsGetType () {
	int pT;

	do {
		pT = 0;
		std::cout << "> Choose form factor type: " << std::endl;
		std::cout << " 1 - Nucleons: p, n " << std::endl;
		std::cout << " 2 - Lambda " << std::endl;
		std::cout << " 3 - Sigma+, Sigma- " << std::endl;
		std::cout << " 4 - Sigma0 " << std::endl;
		std::cout << " 5 - Xi-, Xi0 " << std::endl;
		std::cin >> pT;

		if ( !std::cin ) {
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<int>::max(), '\n');
			std::cout << ">> Error: Not a numeric value!\n" << std::endl;
		} else	if ( pT < 1 || pT > opt ) {
			std::cout << ">> Error: Choose correct value!\n" << std::endl;
		}

	} while ( pT < 1 || pT > opt );
	
	return pT;
}

int HyperonsCheckRange ( const int& pT ) {
	if ( pT < 1 || pT > opt ) { 
		std::cout << ">> Warning: particleType is out of range! Set correct value!\n" << std::endl;
		return 0;
	} else { return 1; };
}
