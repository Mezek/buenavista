/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Debug file.
 */

/**
 * Some useful functions.
 */

/// Show values of FF in the point

int FFinPoint ( char* p, double x )
{
	FFactor debugFF(5);
    debugFF.LoadParameters(p);
	std::cout << "\n t = " << x << std::endl;
	std::cout << "F1S: " << "\t" << debugFF.ScalarOne(x) << std::endl;
	std::cout << "F1V: " << "\t" << debugFF.VectorOne(x) << std::endl;
	std::cout << "F2S: " << "\t" << debugFF.ScalarTwo(x) << std::endl;
	std::cout << "F2V: " << "\t" << debugFF.VectorTwo(x) << std::endl;
	return 0;
}

using namespace ROOT::Minuit2;

/// Perform debug operations.

void performDebug ( char* p, char* f ) {

	std::cout << "\n> Debug parameters:   '" << p << "'" << std::endl;
	std::cout << "> Debug data output:  '" << f << "'" << std::endl;	
	std::cout << "\nDebug stuff:" << std::endl;

	ofstream os (debugFile);

	os << "t = " << -1. << " return: " << FFinPoint ( p, -1. ) << std::endl;
	os << "t = " << -2. << " return: " << FFinPoint ( p, -2. ) << std::endl;
	os << "t = " << -3. << " return: " << FFinPoint ( p, -3. ) << std::endl;
	
	time_t rawtime;
	time ( &rawtime );
	os << ctime (&rawtime) << "Performed debugging." << std::endl;

	os.close();

	std::cout << "\nEnd of debugging." << std::endl;

}
