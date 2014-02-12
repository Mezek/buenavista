/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Experimental data input.
 */

#include "ExperimentalData.h"

namespace ROOT {

	namespace Minuit2 {

/// Read form factor data
		
void ExperimentalData::ReadData ( char* d )
{
	ifstream myDataFile (d);
	//myDataFile.open(dataFile);
	if (myDataFile.is_open()) {
		num = 0;
		double dA,dB,dC;
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
				myDataFile >> dA >> dB >> dC;
				fType.push_back(typeInt);
				fX.push_back(dA);
				fVal.push_back(dB);
				fErrUp.push_back(dC);
				fErrDown.push_back(dC);
				getline (myDataFile,line);
				++num;
			}
		}
		myDataFile.close();
	}
	else std::cerr << "> ReadData: Error! Unable to open data file '" << d << "'!" << std::endl;
}

/// Read form factor data with specific type

void ExperimentalData::ReadData ( char* d, int t )
{
	ifstream myDataFile (d);
	if (myDataFile.is_open()) {
		num = 0;
		double dA,dB,dC;
		while (myDataFile.peek() != EOF) {
			firstChar = myDataFile.peek();
			if ( (firstChar == '%') || (firstChar == '#') || (firstChar == '*') || (firstChar == '/') ) {
				getline (myDataFile,line);
				if ((firstChar == '#')) {
					// Extract data type
					std::string lineSub = line.substr(7);
					// Convert string to int
					std::istringstream stream (lineSub);
					stream >> typeInt;
				}
			}
			else {
				myDataFile >> dA >> dB >> dC;
				if ( typeInt == t) {
					fType.push_back(typeInt);
					fX.push_back(dA);
					fVal.push_back(dB);
					fErrUp.push_back(dC);
					++num;
				}
				getline (myDataFile,line);
			}
		}
		myDataFile.close();
	}
	else std::cerr << "> ReadData: Error! Unable to open data file '" << d << "'!" << std::endl;
	this->typeN(t);
}

/// Read covariance matrix

void ExperimentalData::ReadCovMatrix ( char* d, int n )
{
	ifstream myDataFile (d);
	if (myDataFile.is_open()) {
		num = 0;
		double dA;
		fCov.ResizeTo(n,n);
		while (myDataFile.peek() != EOF) {
			firstChar = myDataFile.peek();
			if ( (firstChar == '%') || (firstChar == '#') || (firstChar == '*') || (firstChar == '/') ) {
				getline (myDataFile,line);
			}
			else {
				for (int i = 0; i < n; ++i) {
					myDataFile >> dA;
					fCov(num,i) = dA;
				}
				getline (myDataFile,line);
				++num;
			}
		}
		myDataFile.close();
	}
	else std::cerr << "> ReadCovMatrix: Error! Unable to open data file '" << d << "'!" << std::endl;
}


/// Show data size

int ExperimentalData::size ()
{
	return fX.size(); 
}

/// Check if type N exists

int ExperimentalData::typeN ( int n )
{
	int k = 0;
	for (int i = 0; i < this->size(); ++i) {
		if ( fType[i] == n ) { ++k; }
	}
	//if ( k == 0 ) { std::cout << "> Warning: Data type '" << n << "' doesn't exist!" << std::endl; }
	return k; 
}

/// Remove specific data and resize data

void ExperimentalData::RemoveDataType ( int n )
{
	int oldSize = this->size();
	int newSize = oldSize - this->typeN(n);
	//std::cout << "> Data: " << oldSize << " " << newSize << std::endl;
	if (oldSize != newSize) {
		int j = 0;
		for (int i = 0; i < oldSize; ++i) {
			if (fType[i] != n) {
				fType[j] = fType[i];
				fX[j] = fX[i];
				fVal[j] = fVal[i];
				fErrUp[j] = fErrUp[i];
				fErrDown[j] = fErrDown[i];
			}
			++j;
		}
		fType.resize(newSize);
		fX.resize(newSize);
		fVal.resize(newSize);
		fErrUp.resize(newSize);
		fErrDown.resize(newSize);
		std::cout << "> Data type: '" << n << "' removed." << std::endl;
	} else {
		std::cout << "> RemoveDataType: Warning! Not such data type: '" << n << "'." << std::endl;
	}
}

/// Check asymmetric errors

void ExperimentalData::CheckData ()
{
	double delta, dVal, dm;
	int k = 0;
	int l = this->size();

	std::cout << "> Data check:" << std::endl;
	std::cout << "Number of checked points:      " << l << std::endl;
	
	for (int i = 0; i < l; ++i) {
		dm = fErrUp[i]-fErrDown[i];
		if ( dm != 0. ) { ++k; } 
		delta = (fErrUp[i] + fErrDown[i])/2.;
		dVal = fVal[i] - delta;
		if ( dVal < 0. ) {
			std::cout << "> CheckData: Warning! Type=" << fType[i] << " data at " << fX[i] 
			<< " : lower error < zero: " << dVal << std::endl;
		}
	}

	std::cout << "Asymmetric errors:             " << k << "\t(" << floor(100.*k/l+0.5) << "%)" << std::endl;
}

/// Show data info

void ExperimentalData::DataInfo ()
{

	int k[140];
	int numCS = 0;
	int numAll = 0;
	
	for (int n = 1; n < 140; ++n) {
		k[n] = this->typeN(n);
		numAll = numAll+k[n];
		if ( n >= 1 ) { numCS = numCS+k[n]; }
	}

	std::cout << "\n> Experimental data:" << std::endl;
	std::cout << "> Data type > 0:                 " << numCS << std::endl;
	std::cout << "> Number of available data:      " << numAll << std::endl;

	if ( this->size() != numAll ) {
		std::cout << "> DataInfo: Warning! Check number of data points!" << std::endl;
	}
}

/// Print all data

void ExperimentalData::ShowData ()
{

	int l = this->size();	

	std::cout << "\n> Experimental data:" << std::endl;
	std::cout << "  No.     Type        X    Value  ErrorUp ErrorDown" << std::endl;

	for (int i = 0; i < l; ++i) {
		std::cout.width(4);
		std::cout << i+1 << ". ";
		std::cout.width(8);
		std::cout << fType[i] << " ";
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
