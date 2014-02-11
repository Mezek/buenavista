/**
 * $Date$
 * $Revision$
 * $Author$
 * $Header$
 * $Id$
 *
 * @file
 * @brief	Header for experimental data input.
 */

#ifndef _ExperimentalData_H_
#define _ExperimentalData_H_

#include "TMatrixD.h"

namespace ROOT {

	namespace Minuit2 {
		
class ExperimentalData {

  private:
		
	char firstChar;
	int typeInt;
	std::string line;
	int num;
	std::vector<int> fType;
	std::vector<double> fX;
	std::vector<double> fVal;
	std::vector<double> fErrUp;
	std::vector<double> fErrDown;
	std::vector<double> fTheta;
	std::vector<double> fEnergy;
	TMatrixD fCov;
	
  public:
		
	std::vector<int> Type() const { return fType; }
	std::vector<double> X() const { return fX; }
	std::vector<double> Val() const { return fVal; }
	std::vector<double> ErrUp() const { return fErrUp; }
	std::vector<double> ErrDown() const { return fErrDown; }
	std::vector<double> Theta() const { return fTheta; }
	std::vector<double> Energy() const { return fEnergy; }
	TMatrixD Cov() const { return fCov; }

	ExperimentalData () { typeInt = 0; }; /**< default type=0 */
	void ReadData ( char* );
	void ReadData ( char*, int );
	void ReadDataCrossSection ( char* );
	void ReadCovMatrix ( char*, int );
	void RemoveDataType ( int );
	void CheckData ();
	void DataInfo ();
	void ShowData ();
	int size ();
	int typeN ( int );
	~ExperimentalData () {};
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // _ExperimentalData_H_
