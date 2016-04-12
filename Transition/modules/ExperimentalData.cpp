/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/ExperimentalData.cpp $
 * $Id$
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
					sSeries.push_back(typeInt);
				}
				if ((firstChar == '*')) {
					// Extract data type
					std::string lineSub = line.substr(7);
					//std::cout << ">> " << lineSub << " " << std::endl;
					sName.push_back(lineSub);
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

/// Read form factor data with specific type

void ExperimentalData::ReadData ( char* d, int t )
{
	ifstream myDataFile (d);
	if (myDataFile.is_open()) {
		num = 0;
		double dA,dB,dC,dD;
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
				myDataFile >> dA >> dB >> dC >> dD;
				if ( typeInt == t) {
					fType.push_back(typeInt);
					fX.push_back(dA);
					fVal.push_back(dB);
					fErrUp.push_back(dC);
					fErrDown.push_back(dD);
					fTheta.push_back(0.);
					fEnergy.push_back(0.);
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

/// Read cross section data

void ExperimentalData::ReadDataCrossSection ( char* d )
{
	ifstream myDataFile (d);
	if (myDataFile.is_open()) {
		num = 0;
		double dA,dB,dC,dD,dE;
		while (myDataFile.peek() != EOF) {
			firstChar = myDataFile.peek();
			if ( (firstChar == '%') || (firstChar == '#') || (firstChar == '*') || (firstChar == '/') ) {
				getline (myDataFile,line);
				if ((firstChar == '#')) {
					std::string lineSub = line.substr(7);
					std::istringstream stream (lineSub);
					stream >> typeInt;
				}
			}
			else {
				myDataFile >> dA >> dB >> dC >> dD >> dE;
				fType.push_back(typeInt);
				fX.push_back(dA);
				fVal.push_back(dB);
				fErrUp.push_back(dC);
				fErrDown.push_back(dC);
				fTheta.push_back(dD);
				fEnergy.push_back(dE);
				getline (myDataFile,line);
				++num;
			}
		}
		myDataFile.close();
	}
	else std::cerr << "> ReadDataCrossSection: Error! Unable to open data file '" << d << "'!" << std::endl;
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
				fTheta[j] = fTheta[i];
				fEnergy[j] = fEnergy[i];
			}
			++j;
		}
		fType.resize(newSize);
		fX.resize(newSize);
		fVal.resize(newSize);
		fErrUp.resize(newSize);
		fErrDown.resize(newSize);
		fTheta.resize(newSize);
		fEnergy.resize(newSize);
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
		if ( fType[i] == 7 ) {
			dVal = -fVal[i] - delta;
		} else {
			dVal = fVal[i] - delta;
		}
		if ( dVal < 0. ) {
			std::cout << "> CheckData: Warning! Type=";
			std::cout << fType[i] << " data at ";
			std::cout.width(7);
			std::cout << fX[i]<< " : y- error exceed 0.0: " << dVal << std::endl;
		}
	}

	std::cout << "Asymmetric errors:             " << k << "\t(" << std::floor(100.*k/l+0.5) << "%)" << std::endl;
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
		if ( n >= 100 ) { numCS = numCS+k[n]; }
	}

	std::cout << "\n> Experimental data:" << std::endl;
	std::cout << "> Data protonElectric:           " << k[1]+k[2] << std::endl;
	std::cout << "> Data protonMagnetic:           " << k[3]+k[4] << std::endl;
	std::cout << "> Data neutronElectric:          " << k[5]+k[6] << std::endl;
	std::cout << "> Data neutronMagnetic:          " << k[7]+k[8] << std::endl;
	std::cout << "> Data protonRatios:             " << k[9] << std::endl;
	std::cout << "> Data neutronRatios:            " << k[10] << std::endl;
	std::cout << "> Data cross section ratios:     " << numCS << std::endl;
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
	std::cout << "  No.     Type        X    Value  ErrorUp ErrorDown   Theta   Energy" << std::endl;

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
		std::cout << fErrDown[i] << " ";
		std::cout.width(8);
		std::cout << fTheta[i] << " ";
		std::cout.width(8);
		std::cout << fEnergy[i] << std::endl;
	}
}


	}  // namespace Minuit2

}  // namespace ROOT



