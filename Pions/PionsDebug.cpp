/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Debug file.
 */

using namespace ROOT::Minuit2;

/// Perform debug operations.

void performDebug ( char* p, char* f ) {

	std::cout << "\n> Debug parameters:   '" << p << "'" << std::endl;
	std::cout << "> Debug data output:  '" << f << "'" << std::endl;	
	std::cout << "\nDebug stuff:" << std::endl;


	std::cout << "> Values of data: " << std::endl;
	char dataFile[] = "../Data/dataPi3Andrej.dat";	
	ExperimentalData A;
	A.ReadData(dataFile);
	A.DataInfo();

}