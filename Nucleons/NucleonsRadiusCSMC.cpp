/**
 * $Date$
 * $Revision$
 * $Author$
 * 
 * @file
 * @brief	Calculate proton radius value from file, v. MC.
 */

using namespace ROOT::Minuit2;

/// Read radius values from file and calculate final value and uncertainty.

void performRadius ( char* f) {

	std::cout << "\n> Radius data:        " << std::endl;
	std::cout << "> Input data file:     '" << f << "'" << std::endl;

	ifstream inps (f);
	std::vector<double> rval;
	std::string line;
	int num;
	if (inps.is_open()) {
		num = 0;
		double dA;
		while (inps.peek() != EOF) {
			inps >> dA;
			rval.push_back(dA);
			getline (inps,line);
			//std::cout << rval[num] << std::endl;
			++num;
		}
		inps.close();
	}
	else std::cerr << "> Error: Unable to open input file '" << f << "'!" << std::endl;
	

	double sigma = 0.;
	double variance = 0.;
	double mean = 0.;

	for (int  i = 0; i < num; ++i) {
		mean += rval[i];
		//std::cout << mean << std::endl;
	}
	mean = mean/num;

	for (int  i = 0; i < num; ++i) {
		variance += (rval[i]-mean)*(rval[i]-mean);
		//std::cout << variance << std::endl;
	}
	sigma = sqrt(variance/(num-1));

	std::cout << "\nProton radius value : Proton radius uncertainty" << std::endl;
	std::cout.width(17); std::cout << mean;
	std::cout << "  +/- ";
	std::cout.width(17); std::cout << sqrt(sigma) << std::endl;

	std::cout << "\nNumber of samples   : " << num << std::endl;
}
