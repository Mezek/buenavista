/**
 * $Date: $
 * $Revision: $
 * $Author: $
 * $HeadURL: $
 * $Id: $
 *
 * @file
 * @brief	Experimental data input.
 */

#include "DeltaData.h"

namespace ROOT {

	namespace Minuit2 {

/// Read form factor data
		
void DeltaData::ReadData ( const char* d )
{
	ifstream myDataFile (d);
	//myDataFile.open(dataFile);
	if (myDataFile.is_open()) {
		num = 0;
		double dA,dB,dC,dD;
		while (myDataFile.peek() != EOF) {
			firstChar = myDataFile.peek();
			if ( (firstChar == '%') || (firstChar == '#') || (firstChar == '*') || (firstChar == '/') ) {
				getline (myDataFile,line);
				//std::cout << line << std::endl;
				if ((firstChar == '#')) {
					// Extract data type
					std::string lineSub = line.substr(7);
					// Convert string to int
					std::istringstream stream (lineSub);
					stream >> typeInt;
					//std::cout << ">> " << lineSub << " " << typeInt << std::endl;
				}
			}
			else {
				myDataFile >> dA >> dB >> dC >> dD;
				fType.push_back(typeInt);
				fX.push_back(dA);
				fVal.push_back(dB);
				fErrUp.push_back(dC);
				fErrDown.push_back(dD);
				fTheta.push_back(0.);
				fEnergy.push_back(0.);
				//std::cout << ">>> " << num << " " << type[num] << " " << x[num] << " " << val[num] << std::endl;
				getline (myDataFile,line);
				++num;
			}
		}
		myDataFile.close();
	}
	else std::cerr << "> ReadData: Error! Unable to open data file '" << d << "'!" << std::endl;
}

/// Show data size

int DeltaData::size ()
{
	return fX.size(); 
}

/// Print all data

void DeltaData::ShowData ()
{
	std::cout << "\n> Experimental data:" << std::endl;
	std::cout << "  No.      X      Value  ErrorUp ErrorDown" << std::endl;

	for (int i = 0; i < this->size(); ++i) {
		std::cout.width(4);
		std::cout << i+1 << ". ";
		std::cout.width(8);
		std::cout << fX[i] << " ";
		std::cout.width(8);
		std::cout << fVal[i] << " ";
		std::cout.width(8);
		std::cout << fErrUp[i] << " ";
		std::cout.width(8);
		std::cout << fErrDown[i] << std::endl;
	}
}


	}  // namespace Minuit2

}  // namespace ROOT



