/**
 * $Date: 2013-05-29 11:53:18 +0200 (Wed, 29 May 2013) $
 * $Revision: 360 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/modules/ExperimentalData.h $
 * $Id: ExperimentalData.h 360 2013-05-29 09:53:18Z bartos $
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
