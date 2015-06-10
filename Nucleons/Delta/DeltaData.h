/**
 * $Date: $
 * $Revision: $
 * $Author: $
 * $HeadURL: $
 * $Id: $
 *
 * @file
 * @brief	Header for experimental data input.
 */

#ifndef _DeltaData_H_
#define _DeltaData_H_


namespace ROOT {

	namespace Minuit2 {
		
class DeltaData {

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
	
  public:
		
	std::vector<int> Type() const { return fType; }
	std::vector<double> X() const { return fX; }
	std::vector<double> Val() const { return fVal; }
	std::vector<double> ErrUp() const { return fErrUp; }
	std::vector<double> ErrDown() const { return fErrDown; }
	std::vector<double> Theta() const { return fTheta; }
	std::vector<double> Energy() const { return fEnergy; }

	DeltaData () {};
	void ReadData ( const char* );
	void ShowData ();
	int size ();
	~DeltaData () {};
};

	}  // namespace Minuit2

}  // namespace ROOT

#endif // _DeltaData_H_
